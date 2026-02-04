import Foundation
import SwiftUI

@MainActor
final class AppModel: ObservableObject {
    @Published var availableModules: [ModuleSpec] = []

    // Pipeline canvas state
    @Published var nodes: [PipelineNode] = []
    @Published var links: [PipelineLink] = []
    @Published var selectedNodeId: UUID?

    // Data + settings
    @Published var samples: [SampleMetadata] = []
    @Published var outputDirectory: String = ""
    @Published var projectName: String = "scanwr-project"
    @Published var recentProjects: [String] = []

    // Logs
    @Published var logs: [String] = []

    // Progress
    @Published var isRunning: Bool = false
    @Published var progressPercent: Double = 0
    @Published var progressMessage: String = ""

    private let rpc = PythonRPCClient()

    var hasProject: Bool {
        !outputDirectory.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
            && !projectName.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
    }

    var projectPath: URL? {
        let base = outputDirectory.trimmingCharacters(in: .whitespacesAndNewlines)
        let name = projectName.trimmingCharacters(in: .whitespacesAndNewlines)
        if base.isEmpty || name.isEmpty { return nil }
        return URL(fileURLWithPath: base).appendingPathComponent(name)
    }

    func loadRecents() {
        recentProjects = UserDefaults.standard.array(forKey: "scanwr.recentProjects") as? [String] ?? []
    }

    func addRecent(_ path: String) {
        var items = recentProjects
        items.removeAll(where: { $0 == path })
        items.insert(path, at: 0)
        recentProjects = Array(items.prefix(12))
        UserDefaults.standard.set(recentProjects, forKey: "scanwr.recentProjects")
    }

    func removeRecent(_ path: String) {
        recentProjects.removeAll(where: { $0 == path })
        UserDefaults.standard.set(recentProjects, forKey: "scanwr.recentProjects")
    }

    func openProject(at url: URL) {
        let fm = FileManager.default
        var projectURL = url
        if url.lastPathComponent == ".scanwr" {
            projectURL = url.deletingLastPathComponent()
        } else if fm.fileExists(atPath: url.appendingPathComponent(".scanwr").path) {
            projectURL = url
        } else {
            appendLog("ERROR: Selected folder is not a scanwr project (missing .scanwr/).")
            return
        }

        outputDirectory = projectURL.deletingLastPathComponent().path
        projectName = projectURL.lastPathComponent
        addRecent(projectURL.path)

        // Load metadata.txt if present
        let metaURL = projectURL.appendingPathComponent(".scanwr/metadata.txt")
        if let text = try? String(contentsOf: metaURL) {
            samples = Self.parseMetadata(text)
            appendLog("Loaded project: \(projectURL.path)")
        } else {
            samples = []
            appendLog("Opened project (no metadata): \(projectURL.path)")
        }
    }

    func createProject(baseDir: URL, name: String) {
        let clean = name.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !clean.isEmpty else { return }

        let projectURL = baseDir.appendingPathComponent(clean)
        let scanwrURL = projectURL.appendingPathComponent(".scanwr")
        let checkpointsURL = scanwrURL.appendingPathComponent("checkpoints")
        let historyURL = scanwrURL.appendingPathComponent("history")

        do {
            try FileManager.default.createDirectory(at: checkpointsURL, withIntermediateDirectories: true)
            try FileManager.default.createDirectory(at: historyURL, withIntermediateDirectories: true)
            let metaURL = scanwrURL.appendingPathComponent("metadata.txt")
            if !FileManager.default.fileExists(atPath: metaURL.path) {
                try "sample\tgroup\tpath\treader\n".write(to: metaURL, atomically: true, encoding: .utf8)
            }
            outputDirectory = baseDir.path
            projectName = clean
            samples = []
            addRecent(projectURL.path)
            appendLog("Created project: \(projectURL.path)")
        } catch {
            appendLog("ERROR: Create project failed: \(error)")
        }
    }

    private static func parseMetadata(_ text: String) -> [SampleMetadata] {
        var out: [SampleMetadata] = []
        let lines = text.split(whereSeparator: \.isNewline)
        guard !lines.isEmpty else { return [] }
        for (idx, raw) in lines.enumerated() {
            if idx == 0 { continue } // header
            let parts = raw.split(separator: "\t", omittingEmptySubsequences: false).map(String.init)
            if parts.count < 3 { continue }
            let sample = parts[0].trimmingCharacters(in: .whitespacesAndNewlines)
            let group = parts[1].trimmingCharacters(in: .whitespacesAndNewlines)
            let path = parts[2].trimmingCharacters(in: .whitespacesAndNewlines)
            let reader = (parts.count >= 4 ? parts[3] : "").trimmingCharacters(in: .whitespacesAndNewlines)
            if sample.isEmpty || group.isEmpty || path.isEmpty { continue }
            out.append(SampleMetadata(sample: sample, group: group, path: path, reader: reader))
        }
        return out
    }

    func appendLog(_ line: String) {
        logs.append(line)
    }

    func ensureBackendStarted() async {
        if rpc.isRunning { return }
        do {
            try await rpc.start(
                onLog: { [weak self] msg in
                    await self?.appendLog(msg)
                },
                onProgress: { [weak self] ev in
                    await self?.setProgress(ev)
                }
            )
            appendLog("Backend: started")
        } catch {
            appendLog("Backend ERROR: \(error)")
        }
    }

    func setProgress(_ ev: PythonRPCClient.ProgressEvent) {
        progressPercent = max(0, min(1, ev.percent))
        progressMessage = ev.message
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

    func detectReader(for path: String) async -> ReaderSuggestion? {
        await ensureBackendStarted()
        struct Params: Codable { var path: String }
        do {
            let res: ReaderSuggestion = try await rpc.call(method: "detect_reader", params: Params(path: path))
            return res
        } catch {
            appendLog("detect_reader ERROR: \(error)")
            return nil
        }
    }

    func addNode(spec: ModuleSpec, at position: CGPoint) {
        let defaults = defaultParams(for: spec.id)
        nodes.append(PipelineNode(specId: spec.id, position: CGPointCodable(position), params: defaults))
    }

    func defaultParams(for specId: String) -> [String: JSONValue] {
        switch specId {
        case "scanpy.pp.filter_cells":
            return [:]
        case "scanpy.pp.filter_genes":
            return [:]
        default:
            return [:]
        }
    }

    func spec(for specId: String) -> ModuleSpec? {
        availableModules.first(where: { $0.id == specId })
    }

    func nodeBinding(id: UUID) -> Binding<PipelineNode>? {
        guard let idx = nodes.firstIndex(where: { $0.id == id }) else { return nil }
        return Binding(get: { self.nodes[idx] }, set: { self.nodes[idx] = $0 })
    }

    func addLink(from: UUID, to: UUID) {
        guard from != to else { return }
        if links.contains(where: { $0.fromNodeId == from && $0.toNodeId == to }) { return }
        links.append(PipelineLink(fromNodeId: from, toNodeId: to))
    }

    func removeSelectedNode() {
        guard let id = selectedNodeId else { return }
        nodes.removeAll(where: { $0.id == id })
        links.removeAll(where: { $0.fromNodeId == id || $0.toNodeId == id })
        selectedNodeId = nil
    }

    func runPipeline() async {
        let outDir = outputDirectory.trimmingCharacters(in: .whitespacesAndNewlines)
        let projName = projectName.trimmingCharacters(in: .whitespacesAndNewlines)

        guard !samples.isEmpty else { appendLog("ERROR: Add at least one sample first"); return }
        guard !outDir.isEmpty else { appendLog("ERROR: Select an output directory first"); return }
        guard !projName.isEmpty else { appendLog("ERROR: Set a project name"); return }

        let ordered = orderedPipeline()
        guard !ordered.isEmpty else { appendLog("ERROR: Add at least one module"); return }

        await ensureBackendStarted()
        appendLog("Running \(ordered.count) step(s) on \(samples.count) sample(s)…")
        isRunning = true
        progressPercent = 0
        progressMessage = "Starting…"
        defer {
            isRunning = false
        }

        struct RunParams: Codable {
            var outputDir: String
            var projectName: String
            var samples: [SampleMetadata]
            var steps: [PipelineStep]
        }
        struct PipelineStep: Codable {
            var specId: String
            var params: [String: JSONValue]
        }

        do {
            let steps = ordered.map { PipelineStep(specId: $0.specId, params: $0.params) }
            let summary: PipelineRunSummary = try await rpc.call(
                method: "run_pipeline",
                params: RunParams(outputDir: outDir, projectName: projName, samples: samples, steps: steps)
            )
            appendLog("OK: wrote outputs to \(summary.outputDir)")
            for r in summary.results {
                appendLog("OK: \(r.sample) via \(r.reader) final=\(r.finalPath) shape=\(r.shape)")
                for p in r.checkpoints {
                    appendLog("  checkpoint: \(p)")
                }
            }
            progressPercent = 1
            progressMessage = "Done."
        } catch {
            appendLog("run_pipeline ERROR: \(error)")
        }
    }

    // Topological ordering for a simple DAG. For now we enforce a linear chain:
    // - exactly one "start" node (in-degree 0)
    // - each node has at most 1 incoming link and at most 1 outgoing link
    func orderedPipeline() -> [PipelineNode] {
        if nodes.isEmpty { return [] }

        var inCount: [UUID: Int] = Dictionary(uniqueKeysWithValues: nodes.map { ($0.id, 0) })
        var outCount: [UUID: Int] = Dictionary(uniqueKeysWithValues: nodes.map { ($0.id, 0) })
        var next: [UUID: UUID] = [:]

        for l in links {
            inCount[l.toNodeId, default: 0] += 1
            outCount[l.fromNodeId, default: 0] += 1
            // Keep last; we validate counts anyway.
            next[l.fromNodeId] = l.toNodeId
        }

        if inCount.values.filter({ $0 == 0 }).count != 1 { return nodes }
        if inCount.values.contains(where: { $0 > 1 }) { return nodes }
        if outCount.values.contains(where: { $0 > 1 }) { return nodes }

        guard let start = inCount.first(where: { $0.value == 0 })?.key else { return nodes }

        var ordered: [PipelineNode] = []
        var seen: Set<UUID> = []
        var cur: UUID? = start
        while let id = cur, !seen.contains(id) {
            seen.insert(id)
            if let node = nodes.first(where: { $0.id == id }) {
                ordered.append(node)
            }
            cur = next[id]
        }

        if ordered.count != nodes.count { return nodes }
        return ordered
    }
}
