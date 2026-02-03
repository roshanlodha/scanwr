import SwiftUI

struct RootView: View {
    @EnvironmentObject private var model: AppModel
    @State private var showAddModule = false
    @State private var showSettings = false
    @State private var showData = false
    @State private var showConsole = false

    var body: some View {
        VStack(spacing: 0) {
            TopBar(
                onData: { showData = true },
                onAddModule: { showAddModule = true },
                onSettings: { showSettings = true },
                onConsole: { showConsole.toggle() },
                onRun: { Task { await model.runPipeline() } }
            )
            Divider()
            ZStack(alignment: .bottom) {
                CanvasView()
                if showConsole {
                    ConsoleDrawer(lines: model.logs, onClose: { showConsole = false })
                        .transition(.move(edge: .bottom).combined(with: .opacity))
                }
            }
        }
        .task { await model.loadModules() }
        .popover(isPresented: $showData, arrowEdge: .top) {
            DataPopover()
                .frame(width: 860, height: 420)
                .padding(4)
        }
        .sheet(isPresented: $showSettings) { SettingsSheet() }
        .popover(isPresented: $showAddModule, arrowEdge: .top) {
            AddModulePopover { spec in
                // Drop near center; user can drag anywhere.
                model.addNode(spec: spec, at: CGPoint(x: 300, y: 220))
                showAddModule = false
            }
            .frame(width: 420, height: 420)
        }
    }
}

private struct TopBar: View {
    @EnvironmentObject private var model: AppModel

    var onData: () -> Void
    var onAddModule: () -> Void
    var onSettings: () -> Void
    var onConsole: () -> Void
    var onRun: () -> Void

    var body: some View {
        VStack(spacing: 8) {
            HStack(spacing: 10) {
                Button {
                    onData()
                } label: {
                    Label("Data", systemImage: "tray.and.arrow.down")
                }
                .disabled(model.isRunning)

                Button {
                    onAddModule()
                } label: {
                    Label("Add Module", systemImage: "plus.circle")
                }
                .disabled(model.isRunning)

                Spacer()

                Button {
                    onSettings()
                } label: {
                    Image(systemName: "info.circle")
                        .imageScale(.large)
                }
                .help("Settings")
                .disabled(model.isRunning)

                Button {
                    onConsole()
                } label: {
                    Image(systemName: "terminal")
                        .imageScale(.large)
                }
                .help("Console")

                Button {
                    onRun()
                } label: {
                    Image(systemName: "play.circle.fill")
                        .imageScale(.large)
                }
                .help("Run pipeline")
                .disabled(
                    model.isRunning
                    || model.outputDirectory.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
                    || model.projectName.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
                    || model.samples.isEmpty
                    || model.nodes.isEmpty
                )
            }

            if model.isRunning {
                HStack(spacing: 10) {
                    ProgressView(value: model.progressPercent)
                        .frame(maxWidth: .infinity)
                    Text(model.progressMessage)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .lineLimit(1)
                        .frame(width: 260, alignment: .leading)
                }
            }
        }
        .padding(.horizontal, 14)
        .padding(.vertical, 10)
    }
}

private struct ConsoleDrawer: View {
    var lines: [String]
    var onClose: () -> Void

    var body: some View {
        VStack(spacing: 0) {
            HStack {
                Text("Console")
                    .font(.headline)
                Spacer()
                Button {
                    onClose()
                } label: {
                    Image(systemName: "chevron.down")
                        .imageScale(.large)
                }
                .buttonStyle(.plain)
                .help("Close console")
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 10)

            Divider()

            ScrollViewReader { proxy in
                ScrollView {
                    LazyVStack(alignment: .leading, spacing: 4) {
                        ForEach(Array(lines.enumerated()), id: \.offset) { idx, line in
                            Text(line)
                                .font(.system(.caption, design: .monospaced))
                                .foregroundStyle(.primary)
                                .frame(maxWidth: .infinity, alignment: .leading)
                                .id(idx)
                        }
                    }
                    .padding(12)
                    .textSelection(.enabled)
                }
                .background(Color(NSColor.textBackgroundColor))
                .onChange(of: lines.count) { _, _ in
                    if let last = lines.indices.last {
                        proxy.scrollTo(last, anchor: .bottom)
                    }
                }
            }
        }
        .frame(height: 260)
        .background(.ultraThinMaterial)
        .clipShape(RoundedRectangle(cornerRadius: 14))
        .shadow(radius: 14)
        .padding(14)
    }
}

private struct AddModulePopover: View {
    @EnvironmentObject private var model: AppModel
    var onPick: (ModuleSpec) -> Void

    @State private var query: String = ""

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Add a module").font(.headline)
            TextField("Search modules…", text: $query)
            List {
                ForEach(ModuleGroup.allCases, id: \.self) { grp in
                    Section(grp.title) {
                        ForEach(filtered(group: grp)) { spec in
                            Button {
                                onPick(spec)
                            } label: {
                                HStack {
                                    Text(spec.title)
                                    Spacer()
                                    if let ns = spec.namespace, ns != "core" {
                                        Text(ns == "experimental" ? "exp" : "ext")
                                            .font(.caption)
                                            .padding(.horizontal, 6)
                                            .padding(.vertical, 2)
                                            .background(Color.primary.opacity(0.10))
                                            .clipShape(RoundedRectangle(cornerRadius: 6))
                                    }
                                    Text(grp.badge)
                                        .font(.caption)
                                        .padding(.horizontal, 6)
                                        .padding(.vertical, 2)
                                        .background(Color(hex: grp.colorHex).opacity(0.15))
                                        .clipShape(RoundedRectangle(cornerRadius: 6))
                                }
                            }
                        }
                    }
                }
            }
        }
        .padding(12)
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

private extension Color {
    init(hex: String) {
        var s = hex.trimmingCharacters(in: CharacterSet.alphanumerics.inverted)
        if s.count == 3 { s = s.map { "\($0)\($0)" }.joined() }
        var v: UInt64 = 0
        Scanner(string: s).scanHexInt64(&v)
        let r = Double((v >> 16) & 0xFF) / 255.0
        let g = Double((v >> 8) & 0xFF) / 255.0
        let b = Double(v & 0xFF) / 255.0
        self.init(.sRGB, red: r, green: g, blue: b, opacity: 1.0)
    }
}
