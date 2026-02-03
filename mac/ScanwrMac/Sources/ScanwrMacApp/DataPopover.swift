import SwiftUI
import UniformTypeIdentifiers

struct DataPopover: View {
    @EnvironmentObject private var model: AppModel

    @State private var selectedId: UUID?
    @State private var showFolderPicker = false

    @State private var draftSample = ""
    @State private var draftGroup = ""
    @State private var draftPath = ""
    @State private var draftReader = "" // empty means auto
    @State private var detecting = false
    @State private var detectedReason: String = ""

    private let readerOptions: [String] = [
        "", // auto
        "scanpy.read_10x_mtx",
        "scanpy.read_10x_h5",
        "scanpy.read_h5ad",
        "scanpy.read_loom",
        "scanpy.read_mtx",
        "scanpy.read",
    ]

    var body: some View {
        VStack(alignment: .leading, spacing: 12) {
            HStack {
                VStack(alignment: .leading, spacing: 2) {
                    Text("Metadata").font(.headline)
                    Text("Add rows with columns: sample, group, path. Path points to a local folder/file for raw data; Scanpy reader will be auto-detected, but you can override.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
                Spacer()
                Button("Add Row") { addRowFromDraft() }
                    .keyboardShortcut(.defaultAction)
                Button("Remove") { removeSelected() }
                    .disabled(selectedId == nil)
            }

            VStack(alignment: .leading, spacing: 8) {
                HStack(spacing: 10) {
                    TextField("sample", text: $draftSample)
                        .frame(width: 140)
                    TextField("group", text: $draftGroup)
                        .frame(width: 140)

                    HStack(spacing: 6) {
                        TextField("path (folder or file)", text: $draftPath)
                        Button("Browse…") { showFolderPicker = true }
                    }

                    Picker("Reader", selection: $draftReader) {
                        Text("auto").tag("")
                        ForEach(readerOptions.dropFirst(), id: \.self) { r in
                            Text(r).tag(r)
                        }
                    }
                    .frame(width: 190)

                    Button(detecting ? "Detecting…" : "Auto-detect") {
                        Task { await autoDetectReader() }
                    }
                    .disabled(draftPath.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty || detecting)
                }

                if !detectedReason.isEmpty {
                    Text(detectedReason)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
            .padding(10)
            .background(Color.primary.opacity(0.04))
            .clipShape(RoundedRectangle(cornerRadius: 12))

            Table(model.samples, selection: $selectedId) {
                TableColumn("sample", value: \.sample)
                TableColumn("group", value: \.group)
                TableColumn("path") { row in
                    Text(row.path).lineLimit(1)
                }
                TableColumn("reader") { row in
                    Text(row.reader.isEmpty ? "auto" : row.reader)
                        .foregroundStyle(row.reader.isEmpty ? .secondary : .primary)
                }
            }

            HStack {
                Spacer()
                Text("\(model.samples.count) row(s)")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
        }
        .padding(12)
        .fileImporter(
            isPresented: $showFolderPicker,
            allowedContentTypes: [.folder, .data],
            allowsMultipleSelection: false
        ) { result in
            if case .success(let urls) = result, let url = urls.first {
                draftPath = url.path
                Task { await autoDetectReader() }
            }
        }
    }

    private func addRowFromDraft() {
        let s = draftSample.trimmingCharacters(in: .whitespacesAndNewlines)
        let g = draftGroup.trimmingCharacters(in: .whitespacesAndNewlines)
        let p = draftPath.trimmingCharacters(in: .whitespacesAndNewlines)
        if s.isEmpty || g.isEmpty || p.isEmpty { return }

        model.samples.append(SampleMetadata(sample: s, group: g, path: p, reader: draftReader))
        draftSample = ""
        draftGroup = ""
        draftPath = ""
        draftReader = ""
        detectedReason = ""
    }

    private func removeSelected() {
        guard let selectedId else { return }
        model.samples.removeAll(where: { $0.id == selectedId })
        self.selectedId = nil
    }

    private func autoDetectReader() async {
        let p = draftPath.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !p.isEmpty else { return }
        detecting = true
        defer { detecting = false }

        if let res = await model.detectReader(for: p) {
            // Only set reader if user hasn't explicitly chosen something.
            if draftReader.isEmpty {
                draftReader = res.suggested
            }
            detectedReason = "Suggested: \(res.suggested) (\(res.reason))"
        } else {
            detectedReason = "Auto-detect failed."
        }
    }
}
