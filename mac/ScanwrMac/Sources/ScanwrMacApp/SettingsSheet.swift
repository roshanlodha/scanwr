import SwiftUI
import UniformTypeIdentifiers

struct SettingsSheet: View {
    @EnvironmentObject private var model: AppModel
    @Environment(\.dismiss) private var dismiss

    @State private var showingDirPicker = false

    var body: some View {
        VStack(alignment: .leading, spacing: 12) {
            Text("Settings").font(.title3).bold()
            Text("Choose where to save per-step checkpoints and the final output files.")
                .foregroundStyle(.secondary)

            HStack(spacing: 10) {
                TextField("Output directory", text: $model.outputDirectory)
                Button("Browse…") { showingDirPicker = true }
            }

            Spacer()

            HStack {
                Spacer()
                Button("Done") { dismiss() }
                    .keyboardShortcut(.defaultAction)
            }
        }
        .padding(16)
        .frame(width: 640, height: 220)
        .fileImporter(
            isPresented: $showingDirPicker,
            allowedContentTypes: [.folder],
            allowsMultipleSelection: false
        ) { result in
            if case .success(let urls) = result, let url = urls.first {
                model.outputDirectory = url.path
            }
        }
    }
}

