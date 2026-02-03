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
            } else {
                Text("No configuration yet.")
                    .foregroundStyle(.secondary)
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

