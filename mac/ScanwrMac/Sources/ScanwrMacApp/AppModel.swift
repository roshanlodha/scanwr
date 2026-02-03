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

    // Logs
    @Published var logs: [String] = []

    private let rpc = PythonRPCClient()

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
        case "pp.calculate_qc_metrics":
            return [
                "expr_type": .string("counts"),
                "var_type": .string("genes"),
                "qc_vars": .string(""),
                "percent_top": .string("50,100,200,500"),
                "layer": .string(""),
                "use_raw": .bool(false),
                "inplace": .bool(false),
                "log1p": .bool(true),
                "parallel": .string(""),
            ]
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

        guard !samples.isEmpty else { appendLog("ERROR: Add at least one sample first"); return }
        guard !outDir.isEmpty else { appendLog("ERROR: Select an output directory first"); return }

        let ordered = orderedPipeline()
        guard !ordered.isEmpty else { appendLog("ERROR: Add at least one module"); return }

        await ensureBackendStarted()
        appendLog("Running \(ordered.count) step(s) on \(samples.count) sample(s)…")

        struct RunParams: Codable {
            var outputDir: String
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
                params: RunParams(outputDir: outDir, samples: samples, steps: steps)
            )
            appendLog("OK: wrote outputs to \(summary.outputDir)")
            for r in summary.results {
                appendLog("OK: \(r.sample) via \(r.reader) final=\(r.finalPath) shape=\(r.shape)")
                for p in r.checkpoints {
                    appendLog("  checkpoint: \(p)")
                }
            }
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
