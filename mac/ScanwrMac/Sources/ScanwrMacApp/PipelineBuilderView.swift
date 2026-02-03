import SwiftUI

struct PipelineBuilderView: View {
    @EnvironmentObject private var model: AppModel

    @State private var selectedAvailable: ModuleSpec.ID?
    @State private var selectedStep: ModuleInstance.ID?

    var body: some View {
        HStack(spacing: 16) {
            VStack(alignment: .leading, spacing: 10) {
                Text("Modules").font(.title3).bold()
                List(selection: $selectedAvailable) {
                    ForEach(ModuleGroup.allCases, id: \.self) { grp in
                        Section {
                            ForEach(model.availableModules.filter { $0.group == grp }) { spec in
                                HStack {
                                    Text(spec.scanpyQualname)
                                    Spacer()
                                    Text(grp.badge)
                                        .font(.caption)
                                        .padding(.horizontal, 6)
                                        .padding(.vertical, 2)
                                        .background(Color(hex: grp.colorHex).opacity(0.15))
                                        .clipShape(RoundedRectangle(cornerRadius: 6))
                                }
                                .tag(spec.id)
                            }
                        } header: {
                            HStack {
                                Rectangle()
                                    .fill(Color(hex: grp.colorHex))
                                    .frame(width: 4)
                                Text("\(grp.title): \(grp.badge)")
                                    .foregroundStyle(Color(hex: grp.colorHex))
                                    .font(.headline)
                            }
                        }
                    }
                }
                Button("Add to Pipeline") {
                    guard let id = selectedAvailable,
                          let spec = model.availableModules.first(where: { $0.id == id }) else { return }
                    model.addModule(spec)
                }
                .disabled(selectedAvailable == nil)
            }
            .frame(minWidth: 360)

            VStack(alignment: .leading, spacing: 10) {
                Text("Pipeline").font(.title3).bold()
                List(selection: $selectedStep) {
                    ForEach(model.pipeline) { step in
                        let title = model.availableModules.first(where: { $0.id == step.specId })?.scanpyQualname ?? step.specId
                        Text(title)
                            .tag(step.id)
                    }
                    .onDelete { idxs in
                        model.pipeline.remove(atOffsets: idxs)
                    }
                }
                HStack {
                    Button("Move Up") { moveSelected(-1) }
                        .disabled(!canMove(-1))
                    Button("Move Down") { moveSelected(1) }
                        .disabled(!canMove(1))
                    Spacer()
                    Button("Remove") { removeSelected() }
                        .disabled(selectedStep == nil)
                }

                Divider()
                StepSettingsView(selectedStepId: $selectedStep)
            }
            .frame(minWidth: 520)
        }
        .padding(16)
    }

    private func indexOfSelected() -> Int? {
        guard let selectedStep else { return nil }
        return model.pipeline.firstIndex(where: { $0.id == selectedStep })
    }

    private func canMove(_ delta: Int) -> Bool {
        guard let idx = indexOfSelected() else { return false }
        let newIdx = idx + delta
        return newIdx >= 0 && newIdx < model.pipeline.count
    }

    private func moveSelected(_ delta: Int) {
        guard let idx = indexOfSelected() else { return }
        let newIdx = idx + delta
        guard newIdx >= 0 && newIdx < model.pipeline.count else { return }
        let item = model.pipeline.remove(at: idx)
        model.pipeline.insert(item, at: newIdx)
        selectedStep = item.id
    }

    private func removeSelected() {
        guard let idx = indexOfSelected() else { return }
        model.pipeline.remove(at: idx)
        selectedStep = nil
    }
}

private struct StepSettingsView: View {
    @EnvironmentObject private var model: AppModel
    @Binding var selectedStepId: ModuleInstance.ID?

    var body: some View {
        if let selectedStepId,
           let idx = model.pipeline.firstIndex(where: { $0.id == selectedStepId }) {
            let specId = model.pipeline[idx].specId
            if specId == "pp.calculate_qc_metrics" {
                CalculateQCMetricsSettings(step: $model.pipeline[idx])
            } else {
                Text("No settings for \(specId)")
                    .foregroundStyle(.secondary)
            }
        } else {
            Text("Select a step to configure it.")
                .foregroundStyle(.secondary)
        }
    }
}

private struct CalculateQCMetricsSettings: View {
    @Binding var step: ModuleInstance

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Calculate QC Metrics").font(.headline)

            LabeledContent("qc_vars (comma-separated)") {
                TextField("", text: Binding(
                    get: { step.params["qc_vars"]?.stringValue ?? "" },
                    set: { step.params["qc_vars"] = .string($0) }
                ))
                .frame(width: 320)
            }
            LabeledContent("percent_top (e.g. 50,100)") {
                TextField("", text: Binding(
                    get: { step.params["percent_top"]?.stringValue ?? "" },
                    set: { step.params["percent_top"] = .string($0) }
                ))
                .frame(width: 320)
            }
            Toggle("log1p", isOn: Binding(
                get: { step.params["log1p"]?.boolValue ?? true },
                set: { step.params["log1p"] = .bool($0) }
            ))
        }
        .padding(.top, 6)
    }
}

private extension Color {
    init(hex: String) {
        var s = hex.trimmingCharacters(in: CharacterSet.alphanumerics.inverted)
        if s.count == 3 {
            s = s.map { "\($0)\($0)" }.joined()
        }
        var v: UInt64 = 0
        Scanner(string: s).scanHexInt64(&v)
        let r = Double((v >> 16) & 0xFF) / 255.0
        let g = Double((v >> 8) & 0xFF) / 255.0
        let b = Double(v & 0xFF) / 255.0
        self.init(.sRGB, red: r, green: g, blue: b, opacity: 1.0)
    }
}
