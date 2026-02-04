import SwiftUI
import AppKit
import UniformTypeIdentifiers

struct ExploreDataView: View {
    @EnvironmentObject private var model: AppModel
    @Environment(\.dismiss) private var dismiss

    @State private var selectedSample: String = ""
    @State private var inspect: AdataInspectResult?
    @State private var isInspecting: Bool = false
    @State private var inspectError: String?

    @State private var keySource: KeySource = .obs
    @State private var keyFilter: String = ""
    @State private var selectedKeys: [String] = []

    // scanpy.pl.violin params
    @State private var groupby: String = ""
    @State private var log: Bool = false
    @State private var useRaw: TriBool = .defaultValue
    @State private var stripplot: Bool = true
    @State private var jitterMode: JitterMode = .bool
    @State private var jitterBool: Bool = true
    @State private var jitterValue: String = ""
    @State private var size: String = "1"
    @State private var layer: String = ""
    @State private var densityNorm: String = "width"
    @State private var orderCSV: String = ""
    @State private var multiPanel: Bool = false
    @State private var xlabel: String = ""
    @State private var ylabel: String = ""
    @State private var rotation: String = ""
    @State private var scale: String = ""
    @State private var show: TriBool = .defaultValue
    @State private var extraKwds: [String: JSONValue] = [:]

    // Output
    @State private var savePath: String = ""
    @State private var lastSVGURL: URL?
    @State private var plotError: String?
    @State private var isPlotting: Bool = false

    var body: some View {
        VStack(spacing: 0) {
            header
            Divider()
            HSplitView {
                ScrollView(.vertical) {
                    VStack(alignment: .leading, spacing: 14) {
                        dataSection
                        keysSection
                        violinSection
                        outputSection
                    }
                    .padding(12)
                }
                .frame(minWidth: 540, idealWidth: 560)

                previewSection
                    .frame(minWidth: 420, idealWidth: 460)
            }
        }
        .frame(maxWidth: .infinity, maxHeight: .infinity)
        .onAppear {
            if selectedSample.isEmpty, let first = model.samples.first?.sample {
                selectedSample = first
                resetForSample()
            }
            if savePath.isEmpty { setDefaultSavePath() }
        }
    }

    private var header: some View {
        HStack {
            VStack(alignment: .leading, spacing: 2) {
                Text("Explore Data")
                    .font(.title2)
                if let p = model.projectPath?.path {
                    Text(p)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .textSelection(.enabled)
                        .lineLimit(1)
                }
            }
            Spacer()
            Button {
                dismiss()
            } label: {
                Image(systemName: "xmark.circle.fill")
                    .imageScale(.large)
                    .foregroundStyle(.secondary)
            }
            .buttonStyle(.plain)
            .help("Close")
        }
        .padding(12)
    }

    private var dataSection: some View {
        GroupBox("Data") {
            VStack(alignment: .leading, spacing: 10) {
                Picker("Sample", selection: $selectedSample) {
                    Text("Select…").tag("")
                    ForEach(model.samples.map(\.sample), id: \.self) { s in
                        Text(s).tag(s)
                    }
                }
                .onChange(of: selectedSample) { _, _ in
                    resetForSample()
                }

                LabeledContent("h5ad") {
                    Text(h5adPathForSelection() ?? "—")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .textSelection(.enabled)
                        .frame(maxWidth: .infinity, alignment: .leading)
                }

                HStack {
                    Button("Load keys") { Task { await loadKeys() } }
                        .disabled(selectedSample.isEmpty || isInspecting)
                    if isInspecting { ProgressView().controlSize(.small) }
                    Spacer()
                }

                if let inspectError {
                    Text(inspectError)
                        .font(.caption)
                        .foregroundStyle(.red)
                        .textSelection(.enabled)
                }
            }
            .padding(6)
        }
    }

    private var keysSection: some View {
        GroupBox("Keys") {
            VStack(alignment: .leading, spacing: 10) {
                Picker("Source", selection: $keySource) {
                    Text("obs").tag(KeySource.obs)
                    Text("var").tag(KeySource.`var`)
                }
                .pickerStyle(.segmented)
                .disabled(inspect == nil)

                TextField("Filter…", text: $keyFilter)
                    .disabled(inspect == nil)

                HStack(alignment: .top, spacing: 12) {
                    VStack(alignment: .leading, spacing: 6) {
                        Text("Available")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                        KeyPickerList(
                            items: filteredAvailableKeys(),
                            emptyText: "No keys.",
                            onPick: { addKey($0) }
                        )
                        .frame(height: 220)
                    }
                    VStack(alignment: .leading, spacing: 6) {
                        Text("Selected")
                            .font(.caption)
                            .foregroundStyle(.secondary)
                        SelectedKeyList(
                            items: selectedKeys,
                            emptyText: "No keys selected.",
                            onRemove: { k in selectedKeys.removeAll(where: { $0 == k }) }
                        )
                        .frame(height: 220)
                    }
                }

                if let inspect, inspect.varNamesTruncated {
                    Text("Note: showing \(inspect.varNames.count) / \(inspect.varNamesTotal) `var_names` (truncated). Use the filter box to find genes.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
            .padding(6)
        }
    }

    private var violinSection: some View {
        GroupBox("Violin") {
            VStack(alignment: .leading, spacing: 10) {
                if let inspect {
                    Picker("groupby", selection: $groupby) {
                        Text("(none)").tag("")
                        ForEach(inspect.groupbyCandidates, id: \.self) { k in
                            Text(k).tag(k)
                        }
                        if !inspect.groupbyCandidates.contains(groupby), !groupby.isEmpty {
                            Text(groupby).tag(groupby)
                        }
                    }
                } else {
                    TextField("groupby (obs col)", text: $groupby)
                        .disabled(true)
                }

                Toggle("log", isOn: $log)

                Picker("use_raw", selection: $useRaw) {
                    Text("default").tag(TriBool.defaultValue)
                    Text("true").tag(TriBool.trueValue)
                    Text("false").tag(TriBool.falseValue)
                }
                .frame(maxWidth: 240)

                Toggle("stripplot", isOn: $stripplot)

                HStack {
                    Picker("jitter", selection: $jitterMode) {
                        Text("bool").tag(JitterMode.bool)
                        Text("float").tag(JitterMode.float)
                    }
                    .frame(width: 110)
                    if jitterMode == .bool {
                        Toggle("true", isOn: $jitterBool).labelsHidden()
                    } else {
                        TextField("e.g. 0.4", text: $jitterValue)
                            .frame(width: 120)
                    }
                    Spacer()
                }

                TextField("size (int)", text: $size)
                    .frame(maxWidth: 240)

                if let inspect {
                    Picker("layer", selection: $layer) {
                        Text("(default)").tag("")
                        ForEach(inspect.layers, id: \.self) { l in
                            Text(l).tag(l)
                        }
                        if !inspect.layers.contains(layer), !layer.isEmpty {
                            Text(layer).tag(layer)
                        }
                    }
                    .frame(maxWidth: 320)
                } else {
                    TextField("layer", text: $layer)
                        .disabled(true)
                }

                TextField("density_norm (width|area|count)", text: $densityNorm)
                    .frame(maxWidth: 320)

                TextField("order (comma-separated)", text: $orderCSV)

                Toggle("multi_panel", isOn: $multiPanel)

                TextField("xlabel", text: $xlabel)
                TextField("ylabel (or comma-separated)", text: $ylabel)
                TextField("rotation (float)", text: $rotation)
                TextField("scale (optional)", text: $scale)

                Picker("show", selection: $show) {
                    Text("default").tag(TriBool.defaultValue)
                    Text("true").tag(TriBool.trueValue)
                    Text("false").tag(TriBool.falseValue)
                }
                .frame(maxWidth: 240)

                Text("Note: plots are always rendered headlessly; `show` is accepted but ignored for display.")
                    .font(.caption)
                    .foregroundStyle(.secondary)

                KeyValueEditor(values: $extraKwds, title: "kwds")
            }
            .padding(6)
        }
    }

    private var outputSection: some View {
        GroupBox("Output") {
            VStack(alignment: .leading, spacing: 10) {
                LabeledContent("SVG path") {
                    TextField("", text: $savePath)
                        .textSelection(.enabled)
                }

                HStack {
                    Button("Choose…") { chooseSavePath() }
                    Button("Default") { setDefaultSavePath() }
                    Spacer()
                }

                HStack {
                    Button("Plot") { Task { await plot() } }
                        .disabled(selectedKeys.isEmpty || selectedSample.isEmpty || isPlotting)
                    if isPlotting { ProgressView().controlSize(.small) }
                    Spacer()
                }

                if let plotError {
                    Text(plotError)
                        .font(.caption)
                        .foregroundStyle(.red)
                        .textSelection(.enabled)
                }
            }
            .padding(6)
        }
    }

    private var previewSection: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Preview").font(.headline)
            if let url = lastSVGURL {
                SVGWebView(fileURL: url)
                    .frame(maxWidth: .infinity, maxHeight: .infinity)
                    .background(Color(NSColor.textBackgroundColor))
                    .clipShape(RoundedRectangle(cornerRadius: 12))
                    .overlay(
                        RoundedRectangle(cornerRadius: 12)
                            .stroke(Color.secondary.opacity(0.2), lineWidth: 1)
                    )
                HStack {
                    Text(url.path)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .textSelection(.enabled)
                        .lineLimit(1)
                    Spacer()
                    Button("Reveal") { NSWorkspace.shared.activateFileViewerSelecting([url]) }
                }
            } else {
                VStack(spacing: 10) {
                    Spacer()
                    Text("No plot yet.")
                        .foregroundStyle(.secondary)
                    Spacer()
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
            }
        }
        .padding(12)
        .frame(maxWidth: .infinity, maxHeight: .infinity)
    }

    private func resetForSample() {
        inspect = nil
        inspectError = nil
        plotError = nil
        selectedKeys = []
        groupby = ""
        layer = ""
        setDefaultSavePath()
    }

    private func h5adPathForSelection() -> String? {
        guard !selectedSample.isEmpty, let project = model.projectPath else { return nil }
        let safe = AppModel.sanitizeFilename(selectedSample)
        let url = project.appendingPathComponent(".scanwr/checkpoints/\(safe).h5ad")
        return url.path
    }

    private func loadKeys() async {
        plotError = nil
        inspectError = nil
        guard let path = h5adPathForSelection() else { return }
        guard FileManager.default.fileExists(atPath: path) else {
            inspectError = "Missing .h5ad for sample. Run the pipeline first (expected at: \(path))"
            return
        }
        isInspecting = true
        defer { isInspecting = false }
        do {
            let res = try await model.inspectH5ad(path: path, varNamesLimit: 5000)
            inspect = res
            if groupby.isEmpty {
                groupby = res.groupbyCandidates.first ?? ""
            }
        } catch {
            inspectError = String(describing: error)
        }
    }

    private func filteredAvailableKeys() -> [String] {
        guard let inspect else { return [] }
        let base: [String]
        switch keySource {
        case .obs:
            base = inspect.obsColumns
        case .`var`:
            base = inspect.varNames
        }
        let q = keyFilter.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()
        if q.isEmpty { return base.filter { !selectedKeys.contains($0) } }
        return base.filter { $0.lowercased().contains(q) && !selectedKeys.contains($0) }
    }

    private func addKey(_ k: String) {
        if selectedKeys.contains(k) { return }
        selectedKeys.append(k)
    }

    private func setDefaultSavePath() {
        guard let project = model.projectPath else { return }
        guard !selectedSample.isEmpty else { return }
        let safeSample = AppModel.sanitizeFilename(selectedSample)
        let ts = Int(Date().timeIntervalSince1970)
        let url = project.appendingPathComponent(".scanwr/plots/\(safeSample)_violin_\(ts).svg")
        savePath = url.path
    }

    private func chooseSavePath() {
        let panel = NSSavePanel()
        panel.allowedContentTypes = [.svg]
        panel.canCreateDirectories = true
        panel.nameFieldStringValue = savePath.isEmpty ? "violin.svg" : URL(fileURLWithPath: savePath).lastPathComponent
        if let project = model.projectPath {
            panel.directoryURL = project.appendingPathComponent(".scanwr/plots")
        }
        panel.begin { resp in
            if resp == .OK, let url = panel.url {
                savePath = url.path
            }
        }
    }

    private func plot() async {
        plotError = nil
        guard let h5ad = h5adPathForSelection() else { return }
        guard FileManager.default.fileExists(atPath: h5ad) else {
            plotError = "Missing .h5ad for sample. Run the pipeline first."
            return
        }
        if savePath.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            setDefaultSavePath()
        }
        isPlotting = true
        defer { isPlotting = false }

        let jitter: JSONValue? = {
            switch jitterMode {
            case .bool:
                return .bool(jitterBool)
            case .float:
                let d = Double(jitterValue.trimmingCharacters(in: .whitespacesAndNewlines)) ?? 0
                return .number(d)
            }
        }()

        let order: [String]? = {
            let s = orderCSV.trimmingCharacters(in: .whitespacesAndNewlines)
            if s.isEmpty { return nil }
            let parts = s.split(separator: ",").map { $0.trimmingCharacters(in: .whitespacesAndNewlines) }.filter { !$0.isEmpty }
            return parts.isEmpty ? nil : parts
        }()

        let request = ViolinPlotRequest(
            h5adPath: h5ad,
            keys: selectedKeys,
            groupby: groupby.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ? nil : groupby,
            log: log,
            useRaw: useRaw.boolValue,
            stripplot: stripplot,
            jitter: jitter,
            size: Int(size.trimmingCharacters(in: .whitespacesAndNewlines)) ?? 1,
            layer: layer.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ? nil : layer,
            densityNorm: densityNorm.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ? "width" : densityNorm,
            order: order,
            multiPanel: multiPanel,
            xlabel: xlabel,
            ylabel: ylabel.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ? nil : ylabel,
            rotation: Double(rotation.trimmingCharacters(in: .whitespacesAndNewlines)),
            show: show.boolValue,
            scale: scale.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty ? nil : scale,
            kwds: extraKwds,
            outputPath: savePath
        )

        do {
            let res = try await model.plotViolin(req: request)
            let url = URL(fileURLWithPath: res.svgPath)
            lastSVGURL = url
        } catch {
            plotError = String(describing: error)
        }
    }
}

private enum KeySource: Hashable {
    case obs
    case `var`
}

private enum TriBool: Hashable {
    case defaultValue
    case trueValue
    case falseValue

    var boolValue: Bool? {
        switch self {
        case .defaultValue: return nil
        case .trueValue: return true
        case .falseValue: return false
        }
    }
}

private enum JitterMode: Hashable {
    case bool
    case float
}

private struct KeyPickerList: View {
    var items: [String]
    var emptyText: String
    var onPick: (String) -> Void

    var body: some View {
        ScrollView(.vertical) {
            LazyVStack(alignment: .leading, spacing: 6) {
                if items.isEmpty {
                    Text(emptyText)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .padding(.vertical, 8)
                } else {
                    ForEach(items, id: \.self) { k in
                        Button {
                            onPick(k)
                        } label: {
                            Text(k)
                                .frame(maxWidth: .infinity, alignment: .leading)
                                .padding(.vertical, 2)
                        }
                        .buttonStyle(.plain)
                    }
                }
            }
            .padding(8)
        }
        .background(Color.secondary.opacity(0.07))
        .clipShape(RoundedRectangle(cornerRadius: 8))
        .overlay(
            RoundedRectangle(cornerRadius: 8)
                .stroke(Color.secondary.opacity(0.12), lineWidth: 1)
        )
    }
}

private struct SelectedKeyList: View {
    var items: [String]
    var emptyText: String
    var onRemove: (String) -> Void

    var body: some View {
        ScrollView(.vertical) {
            LazyVStack(alignment: .leading, spacing: 6) {
                if items.isEmpty {
                    Text(emptyText)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .frame(maxWidth: .infinity, alignment: .leading)
                        .padding(.vertical, 8)
                } else {
                    ForEach(items, id: \.self) { k in
                        HStack(spacing: 8) {
                            Text(k)
                                .frame(maxWidth: .infinity, alignment: .leading)
                                .padding(.vertical, 2)
                            Button {
                                onRemove(k)
                            } label: {
                                Image(systemName: "xmark.circle.fill")
                                    .foregroundStyle(.secondary)
                            }
                            .buttonStyle(.plain)
                        }
                    }
                }
            }
            .padding(8)
        }
        .background(Color.secondary.opacity(0.07))
        .clipShape(RoundedRectangle(cornerRadius: 8))
        .overlay(
            RoundedRectangle(cornerRadius: 8)
                .stroke(Color.secondary.opacity(0.12), lineWidth: 1)
        )
    }
}

private struct KeyValueEditor: View {
    @Binding var values: [String: JSONValue]
    var title: String

    @State private var newKey: String = ""
    @State private var newValue: String = ""
    @State private var newType: String = "string" // string|number|bool|null

    private let typeOptions = ["string", "number", "bool", "null"]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text(title).font(.subheadline).bold()
            Text("Extra kwargs forwarded to Scanpy (advanced). Empty string values are treated as None.")
                .font(.caption)
                .foregroundStyle(.secondary)

            if values.isEmpty {
                Text("No extra kwargs set.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            } else {
                ForEach(values.keys.sorted(), id: \.self) { k in
                    HStack(spacing: 8) {
                        Text(k)
                            .font(.caption)
                            .frame(width: 140, alignment: .leading)
                        TextField("value", text: Binding(
                            get: { displayValue(values[k] ?? .null) },
                            set: { values[k] = parseValue($0, existing: values[k] ?? .string("")) }
                        ))
                        .font(.caption)
                        Button(role: .destructive) {
                            values.removeValue(forKey: k)
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
            values[key] = .string(newValue)
        case "number":
            if let d = Double(newValue.trimmingCharacters(in: .whitespacesAndNewlines)) {
                values[key] = .number(d)
            } else {
                values[key] = .number(0)
            }
        case "bool":
            let v = newValue.trimmingCharacters(in: .whitespacesAndNewlines).lowercased()
            values[key] = .bool(v == "true" || v == "1" || v == "yes")
        default:
            values[key] = .null
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
