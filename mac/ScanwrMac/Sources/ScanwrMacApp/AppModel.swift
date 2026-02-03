import Foundation

@MainActor
final class AppModel: ObservableObject {
    @Published var samples: [SampleRow] = []
    @Published var availableModules: [ModuleSpec] = []
    @Published var pipeline: [ModuleInstance] = []
    @Published var logs: [String] = []

    private var rpc = PythonRPCClient()

    func appendLog(_ line: String) {
        logs.append(line)
    }

    func ensureBackendStarted() async {
        if rpc.isRunning { return }

        do {
            try await rpc.start { [weak self] msg in
                await self?.appendLog(msg)
            }
            appendLog("Backend: started")
        } catch {
            appendLog("Backend ERROR: \(error)")
        }
    }

    func loadModules() async {
        await ensureBackendStarted()
        do {
            let specs: [ModuleSpec] = try await rpc.call(method: "list_modules", params: [:])
            availableModules = specs
        } catch {
            appendLog("list_modules ERROR: \(error)")
        }
    }

    func addModule(_ spec: ModuleSpec) {
        let defaults: [String: JSONValue]
        switch spec.id {
        case "pp.calculate_qc_metrics":
            defaults = [
                "qc_vars": .string(""),
                "percent_top": .string("50,100"),
                "log1p": .bool(true),
            ]
        default:
            defaults = [:]
        }
        pipeline.append(ModuleInstance(specId: spec.id, params: defaults))
    }

    func runPipeline() async {
        guard !samples.isEmpty else { appendLog("ERROR: No samples"); return }
        guard !pipeline.isEmpty else { appendLog("ERROR: No pipeline steps"); return }

        await ensureBackendStarted()
        appendLog("Running \(pipeline.count) step(s) on \(samples.count) sample(s)…")

        do {
            struct RunParams: Codable {
                var samples: [SampleRow]
                var pipeline: [ModuleInstance]
            }
            struct RunResultRow: Codable {
                var sample: String
                var group: String
                var path: String
                var reader: String
                var shape: [Int]
            }

            let res: [RunResultRow] = try await rpc.call(
                method: "run_pipeline",
                params: RunParams(samples: samples, pipeline: pipeline)
            )
            for r in res {
                appendLog("OK: \(r.sample) (group=\(r.group)) via \(r.reader) shape=\(r.shape)")
            }
            appendLog("Done.")
        } catch {
            appendLog("run_pipeline ERROR: \(error)")
        }
    }
}

