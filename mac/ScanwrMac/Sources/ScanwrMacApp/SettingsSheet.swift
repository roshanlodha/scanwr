import SwiftUI
import UniformTypeIdentifiers

struct SettingsSheet: View {
    @EnvironmentObject private var model: AppModel
    @Environment(\.dismiss) private var dismiss

    @State private var showingDirPicker = false

    var body: some View {
        VStack(alignment: .leading, spacing: 14) {
            Text("Settings").font(.title3).bold()

            GroupBox("Pipeline Builder") {
                VStack(alignment: .leading, spacing: 10) {
                    Toggle("Show Console", isOn: $model.pipelineBuilderShowConsole)

                    HStack {
                        Text("Verbosity")
                        Spacer()
                        Text("\(model.verbosity)")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                            .monospacedDigit()
                    }
                    Slider(
                        value: Binding(
                            get: { Double(model.verbosity) },
                            set: { model.setVerbosity(Int($0.rounded())) }
                        ),
                        in: 0...4,
                        step: 1
                    )
                    Text("Lower = quieter console, higher = more details.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                .padding(10)
            }

            GroupBox("Project") {
                VStack(alignment: .leading, spacing: 10) {
                    Text("Choose a project folder. scGUI writes project state into `.scanwr/` inside this folder.")
                        .foregroundStyle(.secondary)

                    Text("Project name: \(model.projectName)")
                        .font(.caption)
                        .foregroundStyle(.secondary)

                    HStack(spacing: 10) {
                        TextField("Project folder", text: $model.outputDirectory)
                        Button("Browse…") { showingDirPicker = true }
                    }

                    if model.outputDirectory.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
                        Text("Project folder is required.")
                            .font(.caption)
                            .foregroundStyle(.red)
                    }
                }
                .padding(10)
            }

            Spacer()

            HStack {
                Spacer()
                Button("Done") { dismiss() }
                    .keyboardShortcut(.defaultAction)
            }
        }
        .padding(16)
        .frame(width: 680, height: 420)
        .fileImporter(
            isPresented: $showingDirPicker,
            allowedContentTypes: [.folder],
            allowsMultipleSelection: false
        ) { result in
            if case .success(let urls) = result, let url = urls.first {
                model.outputDirectory = url.path
                model.projectName = url.lastPathComponent
            }
        }
    }
}
