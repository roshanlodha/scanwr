import SwiftUI
import UniformTypeIdentifiers

struct WelcomeView: View {
    @EnvironmentObject private var model: AppModel

    @State private var showOpenPicker = false
    @State private var showCreatePicker = false
    @State private var newProjectName: String = ""
    @State private var createBaseDir: URL?

    var body: some View {
        HStack(spacing: 0) {
            // Left: branding
            VStack(alignment: .leading, spacing: 14) {
                Spacer()
                Text("scanwr")
                    .font(.system(size: 44, weight: .bold))
                Text("Single-cell pipelines for non-coders.")
                    .font(.title3)
                    .foregroundStyle(.secondary)
                Spacer()
                Text("v0.0.9+")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            .padding(32)
            .frame(minWidth: 320, maxWidth: 360, maxHeight: .infinity, alignment: .topLeading)
            .background(
                LinearGradient(
                    colors: [Color.accentColor.opacity(0.18), Color.accentColor.opacity(0.05)],
                    startPoint: .topLeading,
                    endPoint: .bottomTrailing
                )
            )

            Divider()

            // Right: actions
            VStack(alignment: .leading, spacing: 16) {
                Text("Welcome")
                    .font(.largeTitle).bold()

                HStack(spacing: 12) {
                    Button {
                        showOpenPicker = true
                    } label: {
                        Label("Open Project…", systemImage: "folder")
                    }

                    Button {
                        showCreatePicker = true
                    } label: {
                        Label("New Project…", systemImage: "plus.rectangle.on.folder")
                    }
                }

                if let base = createBaseDir {
                    VStack(alignment: .leading, spacing: 10) {
                        Text("Create project in:")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                        Text(base.path)
                            .font(.caption)
                            .lineLimit(1)
                            .truncationMode(.middle)

                        HStack(spacing: 10) {
                            TextField("Project name", text: $newProjectName)
                                .frame(width: 280)
                            Button("Create") {
                                model.createProject(baseDir: base, name: newProjectName)
                            }
                            .disabled(newProjectName.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty)
                        }
                    }
                    .padding(12)
                    .background(Color.primary.opacity(0.04))
                    .clipShape(RoundedRectangle(cornerRadius: 12))
                }

                Divider()

                Text("Recent")
                    .font(.headline)

                if model.recentProjects.isEmpty {
                    Text("No recent projects yet.")
                        .foregroundStyle(.secondary)
                } else {
                    List {
                        ForEach(model.recentProjects, id: \.self) { path in
                            HStack {
                                Button {
                                    model.openProject(at: URL(fileURLWithPath: path))
                                } label: {
                                    Text(path)
                                        .lineLimit(1)
                                        .truncationMode(.middle)
                                }
                                .buttonStyle(.plain)

                                Spacer()

                                Button(role: .destructive) {
                                    model.removeRecent(path)
                                } label: {
                                    Image(systemName: "xmark.circle.fill")
                                        .foregroundStyle(.secondary)
                                }
                                .buttonStyle(.plain)
                            }
                        }
                    }
                    .frame(maxHeight: 260)
                }

                Spacer()
            }
            .padding(28)
            .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .topLeading)
        }
        .onAppear { model.loadRecents() }
        .fileImporter(
            isPresented: $showOpenPicker,
            allowedContentTypes: [.folder],
            allowsMultipleSelection: false
        ) { result in
            if case .success(let urls) = result, let url = urls.first {
                model.openProject(at: url)
            }
        }
        .fileImporter(
            isPresented: $showCreatePicker,
            allowedContentTypes: [.folder],
            allowsMultipleSelection: false
        ) { result in
            if case .success(let urls) = result, let url = urls.first {
                createBaseDir = url
                newProjectName = ""
            }
        }
    }
}

