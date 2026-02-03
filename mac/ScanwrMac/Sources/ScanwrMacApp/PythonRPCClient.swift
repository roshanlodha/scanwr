import Foundation

enum RPCError: Error, CustomStringConvertible {
    case backendNotRunning
    case invalidResponse
    case backendError(code: Int, message: String)
    case launchFailed(String)

    var description: String {
        switch self {
        case .backendNotRunning:
            return "Python backend not running"
        case .invalidResponse:
            return "Invalid JSON-RPC response"
        case .backendError(let code, let message):
            return "Backend error (\(code)): \(message)"
        case .launchFailed(let msg):
            return "Launch failed: \(msg)"
        }
    }
}

actor PendingCalls {
    private var nextId: Int = 1
    private var pending: [Int: CheckedContinuation<Data, Error>] = [:]

    func makeId() -> Int {
        let id = nextId
        nextId += 1
        return id
    }

    func add(id: Int, cont: CheckedContinuation<Data, Error>) {
        pending[id] = cont
    }

    func resolve(id: Int, result: Result<Data, Error>) {
        guard let cont = pending.removeValue(forKey: id) else { return }
        switch result {
        case .success(let data):
            cont.resume(returning: data)
        case .failure(let err):
            cont.resume(throwing: err)
        }
    }

    func failAll(_ err: Error) {
        let items = pending
        pending.removeAll()
        for (_, cont) in items {
            cont.resume(throwing: err)
        }
    }
}

final class PythonRPCClient {
    private var process: Process?
    private var stdinHandle: FileHandle?
    private var stdoutHandle: FileHandle?
    private var readTask: Task<Void, Never>?
    private let pending = PendingCalls()
    private let writeQueue = DispatchQueue(label: "scanwr.rpc.write", qos: .userInitiated)
    private var onLog: (@Sendable (String) async -> Void)?

    var isRunning: Bool { process?.isRunning == true }

    func start(onLog: @escaping @Sendable (String) async -> Void) async throws {
        if isRunning { return }
        self.onLog = onLog

        guard let python = Self.resolvePythonExecutable() else {
            throw RPCError.launchFailed("Set SCANWR_PYTHON to a Python with scanpy installed.")
        }
        guard let serverURL = Self.resolveServerScriptURL() else {
            throw RPCError.launchFailed("Missing bundled resource: scanwr_rpc_server.py")
        }

        let inPipe = Pipe()
        let outPipe = Pipe()

        let p = Process()
        p.executableURL = URL(fileURLWithPath: python)
        p.arguments = [serverURL.path]
        p.standardInput = inPipe
        p.standardOutput = outPipe
        p.standardError = outPipe
        p.environment = (ProcessInfo.processInfo.environment).merging(
            ["PYTHONUNBUFFERED": "1"],
            uniquingKeysWith: { _, new in new }
        )

        try p.run()

        process = p
        stdinHandle = inPipe.fileHandleForWriting
        stdoutHandle = outPipe.fileHandleForReading

        startReader(handle: outPipe.fileHandleForReading)

        struct Ping: Codable { let ok: Bool }
        let res: Ping = try await call(method: "ping", params: [:])
        if !res.ok {
            throw RPCError.launchFailed("Backend ping failed")
        }
    }

    func stop() async {
        readTask?.cancel()
        readTask = nil
        process?.terminate()
        process = nil
        stdinHandle = nil
        stdoutHandle = nil
        await pending.failAll(RPCError.backendNotRunning)
    }

    func call<Params: Encodable, Result: Decodable>(method: String, params: Params) async throws -> Result {
        let data = try await callRaw(method: method, params: params)
        return try JSONDecoder().decode(Result.self, from: data)
    }

    func call<Result: Decodable>(method: String, params: [String: String]) async throws -> Result {
        let data = try await callRaw(method: method, params: params)
        return try JSONDecoder().decode(Result.self, from: data)
    }

    private func callRaw<Params: Encodable>(method: String, params: Params) async throws -> Data {
        guard isRunning, let stdinHandle else { throw RPCError.backendNotRunning }

        let id = await pending.makeId()
        let req = RPCRequest(id: id, method: method, params: AnyEncodable(params))
        let payload = try JSONEncoder().encode(req)

        return try await withCheckedThrowingContinuation { cont in
            Task {
                await pending.add(id: id, cont: cont)
                writeQueue.async {
                    stdinHandle.write(payload)
                    stdinHandle.write(Data([0x0A]))
                }
            }
        }
    }

    private func startReader(handle: FileHandle) {
        readTask?.cancel()
        readTask = Task.detached { [weak self] in
            guard let self else { return }
            do {
                for try await line in handle.bytes.lines {
                    await self.handleLine(String(line))
                }
            } catch {
                await self.pending.failAll(error)
            }
        }
    }

    private func handleLine(_ line: String) async {
        guard !line.isEmpty else { return }
        guard let data = line.data(using: .utf8) else { return }

        guard let obj = (try? JSONSerialization.jsonObject(with: data)) as? [String: Any] else {
            await onLog?("backend: \(line)")
            return
        }

        if let method = obj["method"] as? String, obj["id"] == nil {
            if method == "log" {
                let msg = ((obj["params"] as? [String: Any])?["message"] as? String) ?? ""
                await onLog?(msg)
            }
            return
        }

        guard let id = obj["id"] as? Int else { return }

        if let err = obj["error"] as? [String: Any] {
            let code = err["code"] as? Int ?? -1
            let msg = err["message"] as? String ?? "Unknown error"
            await pending.resolve(id: id, result: .failure(RPCError.backendError(code: code, message: msg)))
            return
        }

        if let result = obj["result"] {
            if let resultData = try? JSONSerialization.data(withJSONObject: result) {
                await pending.resolve(id: id, result: .success(resultData))
            } else {
                await pending.resolve(id: id, result: .failure(RPCError.invalidResponse))
            }
        }
    }

    private static func resolvePythonExecutable() -> String? {
        if let env = ProcessInfo.processInfo.environment["SCANWR_PYTHON"], !env.isEmpty {
            return env
        }
        if let bundled = resolveBundledPythonPath(), FileManager.default.isExecutableFile(atPath: bundled) {
            return bundled
        }
        return nil
    }

    private static func resolveServerScriptURL() -> URL? {
        if let url = Bundle.main.url(forResource: "scanwr_rpc_server", withExtension: "py") {
            return url
        }
        if let url = Bundle.module.url(forResource: "scanwr_rpc_server", withExtension: "py") {
            return url
        }
        let exe = URL(fileURLWithPath: CommandLine.arguments.first ?? "")
        let dir = exe.deletingLastPathComponent()
        let candidate = dir.appendingPathComponent("scanwr_rpc_server.py")
        if FileManager.default.fileExists(atPath: candidate.path) {
            return candidate
        }
        return nil
    }

    private static func resolveBundledPythonPath() -> String? {
        if let base = Bundle.main.resourceURL {
            return base.appendingPathComponent("python/bin/python3").path
        }
        if let base = Bundle.module.resourceURL {
            return base.appendingPathComponent("python/bin/python3").path
        }
        return nil
    }
}

private struct RPCRequest: Encodable {
    let id: Int
    let method: String
    let params: AnyEncodable
}

private struct AnyEncodable: Encodable {
    private let encodeFn: (Encoder) throws -> Void
    init<T: Encodable>(_ value: T) {
        self.encodeFn = value.encode
    }
    func encode(to encoder: Encoder) throws { try encodeFn(encoder) }
}

