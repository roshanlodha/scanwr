import SwiftUI

struct NodeInspector: View {
    @EnvironmentObject private var model: AppModel
    @Binding var node: PipelineNode

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack {
                VStack(alignment: .leading, spacing: 2) {
                    Text(model.spec(for: node.specId)?.title ?? node.specId)
                        .font(.headline)
                    Text(model.spec(for: node.specId)?.scanpyQualname ?? node.specId)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                Spacer()
                Button {
                    model.selectedNodeId = nil
                } label: {
                    Image(systemName: "xmark.circle.fill")
                        .foregroundStyle(.secondary)
                }
                .buttonStyle(.plain)
            }

            Divider()

            if node.specId == "pp.calculate_qc_metrics" {
                CalculateQCMetricsInspector(params: $node.params)
            } else if node.specId == "scanpy.pp.calculate_qc_metrics" {
                CalculateQCMetricsInspector(params: $node.params)
            } else if node.specId == "scanpy.pp.filter_cells" {
                FilterCellsInspector(params: $node.params)
            } else if node.specId == "scanpy.pp.filter_genes" {
                FilterGenesInspector(params: $node.params)
            } else {
                GenericParamsInspector(params: $node.params)
            }

            Divider()

            Button(role: .destructive) {
                model.removeSelectedNode()
            } label: {
                Label("Remove module", systemImage: "trash")
            }
        }
        .padding(12)
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

            LabeledContent("expr_type") {
                TextField("", text: bindingString("expr_type"))
                    .frame(width: 220)
            }

            LabeledContent("var_type") {
                TextField("", text: bindingString("var_type"))
                    .frame(width: 220)
            }

            LabeledContent("qc_vars (comma-separated)") {
                TextField("", text: bindingString("qc_vars"))
                    .frame(width: 220)
            }

            LabeledContent("percent_top (e.g. 50,100,200,500)") {
                TextField("", text: bindingString("percent_top"))
                    .frame(width: 220)
            }

            LabeledContent("layer (optional)") {
                TextField("", text: bindingString("layer"))
                    .frame(width: 220)
            }

            Toggle("use_raw", isOn: bindingBool("use_raw", default: false))
            Toggle("inplace", isOn: bindingBool("inplace", default: false))
            Toggle("log1p", isOn: bindingBool("log1p", default: true))

            LabeledContent("parallel (optional)") {
                TextField("", text: bindingString("parallel"))
                    .frame(width: 220)
            }
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
