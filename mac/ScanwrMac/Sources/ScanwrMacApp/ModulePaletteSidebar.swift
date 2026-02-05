import SwiftUI

struct ModulePaletteSidebar: View {
    @EnvironmentObject private var model: AppModel
    @State private var query: String = ""

    var body: some View {
        VStack(spacing: 10) {
            HStack(spacing: 10) {
                Text("Modules")
                    .font(.headline)
                Spacer()
            }

            TextField("Search…", text: $query)
                .textFieldStyle(.roundedBorder)

            modulesBody
                .frame(maxHeight: .infinity)

            GroupBox("Run Settings") {
                VStack(alignment: .leading, spacing: 10) {
                    ToggleRow(
                        title: "Combine Samples",
                        isOn: Binding<Bool>(
                            get: { model.analysisMode == .concat },
                            set: { model.analysisMode = $0 ? .concat : .perSample }
                        )
                    )

                    ToggleRow(title: "Force Rerun", isOn: $model.forceRerun)
                        .help("Ignores cached checkpoints and reruns from step 1.")
                }
                .padding(.top, 2)
            }
            .disabled(model.isRunning)
        }
        .padding(12)
        .frame(maxHeight: .infinity, alignment: .top)
        .panelChrome()
        .padding(10)
    }

    private var modulesBody: some View {
        Group {
            if model.isLoadingModules {
                VStack(spacing: 10) {
                    Spacer()
                    ProgressView("Loading modules…")
                    Spacer()
                }
            } else if model.availableModules.isEmpty {
                VStack(spacing: 10) {
                    Spacer()
                    Text("No modules available.")
                        .font(.headline)
                    if let err = model.moduleLoadError {
                        Text(err)
                            .font(.caption)
                            .foregroundStyle(.secondary)
                            .multilineTextAlignment(.center)
                            .textSelection(.enabled)
                    } else {
                        Text("Check the Console for backend logs and retry.")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                            .multilineTextAlignment(.center)
                    }
                    Button("Retry") { Task { await model.loadModules() } }
                    Spacer()
                }
            } else {
                List {
                    ForEach(ModuleGroup.allCases, id: \.self) { grp in
                        Section(grp.title) {
                            ForEach(filtered(group: grp)) { spec in
                                ModuleRow(spec: spec)
                            }
                        }
                    }
                }
                .listStyle(.inset)
                .scrollContentBackground(.hidden)
            }
        }
    }

    private func filtered(group: ModuleGroup) -> [ModuleSpec] {
        let q = query.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()
        let base = model.availableModules.filter { $0.group == group }
        if q.isEmpty { return base }
        return base.filter { spec in
            spec.title.lowercased().contains(q)
                || spec.scanpyQualname.lowercased().contains(q)
                || (spec.namespace ?? "").lowercased().contains(q)
        }
    }
}

private struct ToggleRow: View {
    var title: String
    @Binding var isOn: Bool

    var body: some View {
        HStack(spacing: 10) {
            Text(title)
                .frame(maxWidth: .infinity, alignment: .leading)
            Toggle("", isOn: $isOn)
                .labelsHidden()
                .toggleStyle(.switch)
        }
        .frame(maxWidth: .infinity)
    }
}

private struct ModuleRow: View {
    @EnvironmentObject private var model: AppModel
    var spec: ModuleSpec

    var body: some View {
        HStack(spacing: 10) {
            VStack(alignment: .leading, spacing: 2) {
                Text(spec.title)
                    .font(.body)
                    .lineLimit(1)
                Text(spec.scanpyQualname)
                    .font(.caption)
                    .foregroundStyle(.secondary)
                    .lineLimit(1)
            }
            Spacer(minLength: 10)
            if spec.group == .custom || spec.namespace == "custom" {
                Text("rl")
                    .font(.caption)
                    .padding(.horizontal, 6)
                    .padding(.vertical, 2)
                    .background(Color(hex: ModuleGroup.custom.colorHex).opacity(0.18))
                    .clipShape(RoundedRectangle(cornerRadius: 6))
            } else if let ns = spec.namespace, ns != "core" {
                Text(ns == "experimental" ? "exp" : "ext")
                    .font(.caption)
                    .padding(.horizontal, 6)
                    .padding(.vertical, 2)
                    .background(Color.primary.opacity(0.10))
                    .clipShape(RoundedRectangle(cornerRadius: 6))
            }
            if spec.group != .custom {
                Text(spec.group.badge)
                    .font(.caption)
                    .padding(.horizontal, 6)
                    .padding(.vertical, 2)
                    .background(Color(hex: spec.group.colorHex).opacity(0.15))
                    .clipShape(RoundedRectangle(cornerRadius: 6))
            }
        }
        .contentShape(Rectangle())
        .listRowBackground(Color.clear)
        .onTapGesture(count: 2) {
            model.appendStep(spec: spec)
        }
        .onDrag {
            NSItemProvider(object: spec.id as NSString)
        }
        .help("Drag into the pipeline to add, or double-click to append")
    }
}
