import SwiftUI

struct ModulePaletteSidebar: View {
    @EnvironmentObject private var model: AppModel
    @Binding var showSettings: Bool
    @Binding var showConsole: Bool
    @State private var query: String = ""

    var body: some View {
        VStack(spacing: 10) {
            HStack(spacing: 10) {
                Text("Modules")
                    .font(.headline)
                Spacer()

                Button {
                    showSettings = true
                } label: {
                    Image(systemName: "info.circle")
                        .imageScale(.medium)
                }
                .buttonStyle(.plain)
                .help("Settings")
                .disabled(model.isRunning)

                Button {
                    showConsole.toggle()
                } label: {
                    Image(systemName: showConsole ? "terminal.fill" : "terminal")
                        .imageScale(.medium)
                }
                .buttonStyle(.plain)
                .help("Console")

                if model.isRunning {
                    Button {
                        Task { await model.stopRun() }
                    } label: {
                        Image(systemName: "stop.circle.fill")
                            .imageScale(.medium)
                    }
                    .buttonStyle(.plain)
                    .help("Stop run")
                } else {
                    Button {
                        model.startRun()
                    } label: {
                        Image(systemName: "play.circle.fill")
                            .imageScale(.medium)
                    }
                    .buttonStyle(.plain)
                    .help("Run pipeline")
                    .disabled(
                        model.outputDirectory.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
                        || model.projectName.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
                        || model.samples.isEmpty
                        || model.nodes.isEmpty
                    )
                }
            }

            TextField("Search…", text: $query)
                .textFieldStyle(.roundedBorder)

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
        .padding(12)
        .frame(maxHeight: .infinity, alignment: .top)
        .panelChrome()
        .padding(10)
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
            if let ns = spec.namespace, ns != "core" {
                Text(ns == "experimental" ? "exp" : "ext")
                    .font(.caption)
                    .padding(.horizontal, 6)
                    .padding(.vertical, 2)
                    .background(Color.primary.opacity(0.10))
                    .clipShape(RoundedRectangle(cornerRadius: 6))
            }
            Text(spec.group.badge)
                .font(.caption)
                .padding(.horizontal, 6)
                .padding(.vertical, 2)
                .background(Color(hex: spec.group.colorHex).opacity(0.15))
                .clipShape(RoundedRectangle(cornerRadius: 6))
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
