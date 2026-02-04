import SwiftUI

struct NodeInspector: View {
    @EnvironmentObject private var model: AppModel
    @Binding var node: PipelineNode

    var body: some View {
        NodeInspectorWrapper(
            title: model.spec(for: node.specId)?.title ?? node.specId,
            subtitle: model.spec(for: node.specId)?.scanpyQualname ?? node.specId,
            onClose: { model.selectedNodeId = nil },
            content: {
                NodeInspectorContent(specId: node.specId, params: $node.params)

                Divider()

                Button(role: .destructive) {
                    model.removeSelectedNode()
                } label: {
                    Label("Remove module", systemImage: "trash")
                }
            }
        )
    }
}

struct NodeInspectorInline: View {
    @EnvironmentObject private var model: AppModel
    @Binding var node: PipelineNode

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            NodeInspectorContent(specId: node.specId, params: $node.params)

            Divider()

            Button(role: .destructive) {
                model.removeNode(id: node.id)
            } label: {
                Label("Remove module", systemImage: "trash")
            }
        }
        .padding(12)
        .background(Color.primary.opacity(0.04))
        .clipShape(RoundedRectangle(cornerRadius: 12))
    }
}

private struct NodeInspectorWrapper<Content: View>: View {
    var title: String
    var subtitle: String
    var onClose: () -> Void
    let content: Content

    init(
        title: String,
        subtitle: String,
        onClose: @escaping () -> Void,
        @ViewBuilder content: () -> Content
    ) {
        self.title = title
        self.subtitle = subtitle
        self.onClose = onClose
        self.content = content()
    }

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack {
                VStack(alignment: .leading, spacing: 2) {
                    Text(title)
                        .font(.headline)
                    Text(subtitle)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                Spacer()
                Button {
                    onClose()
                } label: {
                    Image(systemName: "xmark.circle.fill")
                        .foregroundStyle(.secondary)
                }
                .buttonStyle(.plain)
            }

            Divider()

            content
        }
        .padding(12)
    }
}

private struct NodeInspectorContent: View {
    var specId: String
    @Binding var params: [String: JSONValue]

    var body: some View {
        if specId == "pp.calculate_qc_metrics" {
            CalculateQCMetricsInspector(params: $params)
        } else if specId == "scanpy.pp.calculate_qc_metrics" {
            CalculateQCMetricsInspector(params: $params)
        } else if specId == "rapids_singlecell.pp.calculate_qc_metrics" {
            CalculateQCMetricsInspector(params: $params)
        } else if specId == "scanpy.pp.filter_cells" {
            FilterCellsInspector(params: $params)
        } else if specId == "rapids_singlecell.pp.filter_cells" {
            FilterCellsInspector(params: $params)
        } else if specId == "scanpy.pp.filter_genes" {
            FilterGenesInspector(params: $params)
        } else if specId == "rapids_singlecell.pp.filter_genes" {
            FilterGenesInspector(params: $params)
        } else if specId == "scanpy.pp.scrublet" {
            ScrubletInspector(params: $params)
        } else if specId == "rapids_singlecell.pp.scrublet" {
            ScrubletInspector(params: $params)
        } else if specId == "scanpy.pp.highly_variable_genes" {
            HighlyVariableGenesInspector(params: $params)
        } else if specId == "rapids_singlecell.pp.highly_variable_genes" {
            HighlyVariableGenesInspector(params: $params)
        } else if specId == "scanpy.pp.normalize_total" {
            NormalizeTotalInspector(params: $params)
        } else if specId == "rapids_singlecell.pp.normalize_total" {
            NormalizeTotalInspector(params: $params)
        } else if specId == "scanpy.tl.pca" {
            PCAInspector(params: $params)
        } else if specId == "rapids_singlecell.pp.pca" {
            PCAInspector(params: $params)
        } else if specId == "scanpy.tl.leiden" {
            LeidenInspector(params: $params)
        } else if specId == "rapids_singlecell.tl.leiden" {
            LeidenInspector(params: $params)
        } else if specId == "scanpy.tl.umap" {
            UMAPInspector(params: $params)
        } else if specId == "rapids_singlecell.tl.umap" {
            UMAPInspector(params: $params)
        } else if specId == "scanpy.tl.rank_genes_groups" {
            RankGenesGroupsInspector(params: $params)
        } else if specId == "rapids_singlecell.tl.rank_genes_groups" {
            RankGenesGroupsInspector(params: $params)
        } else {
            GenericParamsInspector(params: $params)
        }
    }
}

private struct GenericParamsInspector: View {
    @Binding var params: [String: JSONValue]

    @State private var newKey: String = ""
    @State private var newValue: String = ""
    @State private var newType: String = "string" // string|number|bool|null

    private let typeOptions = ["string", "number", "bool", "null"]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Many Scanpy functions need extra args. Add kwargs here (advanced). Empty string values are treated as None.")
                .font(.caption)
                .foregroundStyle(.secondary)

            if params.isEmpty {
                Text("No parameters set.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            } else {
                ForEach(params.keys.sorted(), id: \.self) { k in
                    HStack(spacing: 8) {
                        Text(k)
                            .font(.caption)
                            .frame(width: 140, alignment: .leading)
                        TextField("value", text: Binding(
                            get: { displayValue(params[k] ?? .null) },
                            set: { params[k] = parseValue($0, existing: params[k] ?? .string("")) }
                        ))
                        .font(.caption)
                        Button(role: .destructive) {
                            params.removeValue(forKey: k)
                        } label: {
                            Image(systemName: "trash")
                        }
                        .buttonStyle(.borderless)
                    }
                }
            }

            Divider()

            HStack(spacing: 8) {
                TextField("key", text: $newKey)
                    .frame(width: 180)
                Picker("type", selection: $newType) {
                    ForEach(typeOptions, id: \.self) { t in
                        Text(t).tag(t)
                    }
                }
                .frame(width: 110)
                TextField("value", text: $newValue)
                Button("Add") { addNewParam() }
                    .disabled(newKey.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty)
            }
            .font(.caption)
        }
    }

    private func addNewParam() {
        let key = newKey.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !key.isEmpty else { return }
        switch newType {
        case "string":
            params[key] = .string(newValue)
        case "number":
            if let d = Double(newValue.trimmingCharacters(in: .whitespacesAndNewlines)) {
                params[key] = .number(d)
            } else {
                params[key] = .number(0)
            }
        case "bool":
            let v = newValue.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()
            params[key] = .bool(v == "true" || v == "1" || v == "yes")
        default:
            params[key] = .null
        }
        newKey = ""
        newValue = ""
        newType = "string"
    }

    private func displayValue(_ v: JSONValue) -> String {
        switch v {
        case .string(let s): return s
        case .number(let n):
            if n.rounded() == n { return String(Int(n)) }
            return String(n)
        case .bool(let b): return b ? "true" : "false"
        case .null: return ""
        }
    }

    private func parseValue(_ s: String, existing: JSONValue) -> JSONValue {
        switch existing {
        case .number:
            return .number(Double(s.trimmingCharacters(in: .whitespacesAndNewlines)) ?? 0)
        case .bool:
            let v = s.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()
            return .bool(v == "true" || v == "1" || v == "yes")
        case .null:
            return .string(s)
        case .string:
            return .string(s)
        }
    }
}

private struct CalculateQCMetricsInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()

            Text("Select which gene sets to calculate QC metrics for.")
                .font(.caption)
                .foregroundStyle(.secondary)

            Toggle("mitochondrial (mt)", isOn: bindingBool("use_mt", default: true))
            Toggle("ribosomal (ribo)", isOn: bindingBool("use_ribo", default: true))
            Toggle("hemoglobin (hb)", isOn: bindingBool("use_hb", default: true))

            LabeledContent("percent_top (optional; e.g. 50,100)") {
                TextField("", text: bindingString("percent_top"))
                    .frame(width: 220)
            }

            Toggle("log1p", isOn: bindingBool("log1p", default: true))
        }
    }

    private func bindingString(_ key: String) -> Binding<String> {
        Binding<String>(
            get: { params[key]?.stringValue ?? "" },
            set: { params[key] = $0.isEmpty ? .string("") : .string($0) }
        )
    }

    private func bindingBool(_ key: String, default def: Bool) -> Binding<Bool> {
        Binding<Bool>(
            get: { params[key]?.boolValue ?? def },
            set: { params[key] = .bool($0) }
        )
    }
}

private struct ScrubletInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("If set, `batch_key` is an obs column used to run Scrublet per-batch.")
                .font(.caption)
                .foregroundStyle(.secondary)

            LabeledContent("batch_key") {
                TextField("", text: bindingString("batch_key"))
                    .frame(width: 220)
            }
        }
    }

    private func bindingString(_ key: String) -> Binding<String> {
        Binding<String>(
            get: { params[key]?.stringValue ?? "" },
            set: { params[key] = .string($0) }
        )
    }
}

private struct NormalizeTotalInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Leave `target_sum` blank to use Scanpy’s default.")
                .font(.caption)
                .foregroundStyle(.secondary)

            LabeledContent("target_sum") {
                TextField("", text: bindingString("target_sum"))
                    .frame(width: 220)
            }
        }
    }

    private func bindingString(_ key: String) -> Binding<String> {
        Binding<String>(
            get: { params[key]?.stringValue ?? "" },
            set: { params[key] = .string($0) }
        )
    }
}

private struct HighlyVariableGenesInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Leave values blank to use Scanpy’s defaults.")
                .font(.caption)
                .foregroundStyle(.secondary)

            LabeledContent("n_top_genes") {
                TextField("", text: bindingString("n_top_genes"))
                    .frame(width: 220)
            }

            LabeledContent("batch_key") {
                TextField("", text: bindingString("batch_key"))
                    .frame(width: 220)
            }
        }
    }

    private func bindingString(_ key: String) -> Binding<String> {
        Binding<String>(
            get: { params[key]?.stringValue ?? (params[key]?.doubleValue.map { String(Int($0)) } ?? "") },
            set: { params[key] = .string($0) }
        )
    }
}

private struct LeidenInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Resolution parameter for Leiden clustering.")
                .font(.caption)
                .foregroundStyle(.secondary)

            LabeledContent("res") {
                TextField("", text: bindingString("res"))
                    .frame(width: 220)
            }
        }
    }

    private func bindingString(_ key: String) -> Binding<String> {
        Binding<String>(
            get: { params[key]?.stringValue ?? (params[key]?.doubleValue.map { String($0) } ?? "") },
            set: { params[key] = .string($0) }
        )
    }
}

private struct PCAInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Runs PCA (writes embeddings into `adata.obsm['X_pca']`). Optional plot saved as `Project/plots/{sample}/pca_scatter.png`.")
                .font(.caption)
                .foregroundStyle(.secondary)
            Toggle("Create Scatterplot", isOn: bindingBool("create_scatterplot", default: true))
        }
    }

    private func bindingBool(_ key: String, default def: Bool) -> Binding<Bool> {
        Binding<Bool>(
            get: { params[key]?.boolValue ?? def },
            set: { params[key] = .bool($0) }
        )
    }
}

private struct UMAPInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Runs UMAP (writes embeddings into `adata.obsm['X_umap']`). Optional plot saved as `Project/plots/{sample}/umap_scatter.png`.")
                .font(.caption)
                .foregroundStyle(.secondary)
            Toggle("Create Scatterplot", isOn: bindingBool("create_scatterplot", default: true))
        }
    }

    private func bindingBool(_ key: String, default def: Bool) -> Binding<Bool> {
        Binding<Bool>(
            get: { params[key]?.boolValue ?? def },
            set: { params[key] = .bool($0) }
        )
    }
}

private struct RankGenesGroupsInspector: View {
    @Binding var params: [String: JSONValue]

    private enum Groupby: String, CaseIterable, Identifiable {
        case leiden
        case kmeans
        case louvain

        var id: String { rawValue }

        var title: String {
            switch self {
            case .leiden: return "leiden"
            case .kmeans: return "k-means"
            case .louvain: return "louvain"
            }
        }
    }

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Ranks marker genes per cluster/group (writes results into `adata.uns['rank_genes_groups']`). Optional plots are saved under `Project/plots/{sample}/` as PNG.")
                .font(.caption)
                .foregroundStyle(.secondary)

            LabeledContent("groupby") {
                Picker("", selection: bindingGroupby()) {
                    ForEach(Groupby.allCases) { g in
                        Text(g.title).tag(g.rawValue)
                    }
                }
                .frame(width: 220)
            }

            Toggle("Create Dotplot", isOn: bindingBool("create_dotplot", default: false))
            Toggle("Create Heatmap", isOn: bindingBool("create_heatmap", default: true))
        }
    }

    private func bindingGroupby() -> Binding<String> {
        Binding<String>(
            get: {
                let raw = params["groupby"]?.stringValue ?? "leiden"
                let s = raw.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()
                if s == "k-means" || s == "k_means" { return "kmeans" }
                if Groupby.allCases.map(\.rawValue).contains(s) { return s }
                return "leiden"
            },
            set: { params["groupby"] = .string($0) }
        )
    }

    private func bindingBool(_ key: String, default def: Bool) -> Binding<Bool> {
        Binding<Bool>(
            get: { params[key]?.boolValue ?? def },
            set: { params[key] = .bool($0) }
        )
    }
}

private struct FilterCellsInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Only the thresholds are configurable here; `data` is the running AnnData, and `inplace=True`.")
                .font(.caption)
                .foregroundStyle(.secondary)

            LabeledContent("min_counts") { intField("min_counts") }
            LabeledContent("min_genes") { intField("min_genes") }
            LabeledContent("max_counts") { intField("max_counts") }
            LabeledContent("max_genes") { intField("max_genes") }
        }
    }

    private func intField(_ key: String) -> some View {
        TextField("", text: Binding(
            get: { displayInt(params[key]) },
            set: { setOptionalInt(key: key, text: $0) }
        ))
        .frame(width: 180)
    }

    private func displayInt(_ v: JSONValue?) -> String {
        guard let v else { return "" }
        if case .number(let n) = v { return String(Int(n)) }
        if case .string(let s) = v { return s }
        return ""
    }

    private func setOptionalInt(key: String, text: String) {
        let t = text.trimmingCharacters(in: .whitespacesAndNewlines)
        if t.isEmpty {
            params.removeValue(forKey: key)
            return
        }
        if let i = Int(t) {
            params[key] = .number(Double(i))
        } else {
            params[key] = .string(t)
        }
    }
}

private struct FilterGenesInspector: View {
    @Binding var params: [String: JSONValue]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Parameters").font(.subheadline).bold()
            Text("Only the thresholds are configurable here; `data` is the running AnnData, and `inplace=True`.")
                .font(.caption)
                .foregroundStyle(.secondary)

            LabeledContent("min_counts") { intField("min_counts") }
            LabeledContent("min_cells") { intField("min_cells") }
            LabeledContent("max_counts") { intField("max_counts") }
            LabeledContent("max_cells") { intField("max_cells") }
        }
    }

    private func intField(_ key: String) -> some View {
        TextField("", text: Binding(
            get: { displayInt(params[key]) },
            set: { setOptionalInt(key: key, text: $0) }
        ))
        .frame(width: 180)
    }

    private func displayInt(_ v: JSONValue?) -> String {
        guard let v else { return "" }
        if case .number(let n) = v { return String(Int(n)) }
        if case .string(let s) = v { return s }
        return ""
    }

    private func setOptionalInt(key: String, text: String) {
        let t = text.trimmingCharacters(in: .whitespacesAndNewlines)
        if t.isEmpty {
            params.removeValue(forKey: key)
            return
        }
        if let i = Int(t) {
            params[key] = .number(Double(i))
        } else {
            params[key] = .string(t)
        }
    }
}
