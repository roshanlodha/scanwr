import SwiftUI

struct SamplesView: View {
    @EnvironmentObject private var model: AppModel
    @State private var showingAdd = false

    var body: some View {
        VStack(alignment: .leading, spacing: 12) {
            HStack {
                Text("Samples").font(.title2).bold()
                Spacer()
                Button("Add Sample") { showingAdd = true }
            }

            Table(model.samples) {
                TableColumn("sample", value: \.sample)
                TableColumn("group", value: \.group)
                TableColumn("path", value: \.path)
            }
        }
        .padding(16)
        .sheet(isPresented: $showingAdd) {
            AddSampleSheet { row in
                model.samples.append(row)
            }
        }
    }
}

private struct AddSampleSheet: View {
    var onAdd: (SampleRow) -> Void

    @Environment(\.dismiss) private var dismiss
    @State private var sample = ""
    @State private var group = ""
    @State private var path = ""
    @State private var showPicker = false

    var body: some View {
        VStack(alignment: .leading, spacing: 12) {
            Text("Add Sample").font(.title3).bold()

            Grid(alignment: .leading, horizontalSpacing: 12, verticalSpacing: 10) {
                GridRow {
                    Text("sample")
                    TextField("", text: $sample).frame(width: 360)
                }
                GridRow {
                    Text("group")
                    TextField("", text: $group).frame(width: 360)
                }
                GridRow {
                    Text("path")
                    HStack {
                        TextField("", text: $path).frame(width: 320)
                        Button("Browse…") { showPicker = true }
                    }
                }
            }

            HStack {
                Spacer()
                Button("Cancel") { dismiss() }
                Button("Add") {
                    onAdd(SampleRow(sample: sample.trimmingCharacters(in: .whitespacesAndNewlines),
                                   group: group.trimmingCharacters(in: .whitespacesAndNewlines),
                                   path: path.trimmingCharacters(in: .whitespacesAndNewlines)))
                    dismiss()
                }
                .disabled(sample.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ||
                          group.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ||
                          path.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty)
            }
        }
        .padding(16)
        .frame(minWidth: 520)
        .fileImporter(
            isPresented: $showPicker,
            allowedContentTypes: [.folder, .data],
            allowsMultipleSelection: false
        ) { result in
            if case .success(let urls) = result, let url = urls.first {
                path = url.path
            }
        }
    }
}

