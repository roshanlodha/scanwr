import SwiftUI

struct RootView: View {
    @EnvironmentObject private var model: AppModel
    @State private var showAddModule = false
    @State private var showSettings = false
    @State private var showData = false

    var body: some View {
        VStack(spacing: 0) {
            TopBar(
                onData: { showData = true },
                onAddModule: { showAddModule = true },
                onSettings: { showSettings = true },
                onRun: { Task { await model.runPipeline() } }
            )
            Divider()
            CanvasView()
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
    var onData: () -> Void
    var onAddModule: () -> Void
    var onSettings: () -> Void
    var onRun: () -> Void

    var body: some View {
        HStack(spacing: 10) {
            Button {
                onData()
            } label: {
                Label("Data", systemImage: "tray.and.arrow.down")
            }

            Button {
                onAddModule()
            } label: {
                Label("Add Module", systemImage: "plus.circle")
            }

            Spacer()

            Button {
                onSettings()
            } label: {
                Image(systemName: "info.circle")
                    .imageScale(.large)
            }
            .help("Settings")

            Button {
                onRun()
            } label: {
                Image(systemName: "play.circle.fill")
                    .imageScale(.large)
            }
            .help("Run pipeline")
        }
        .padding(.horizontal, 14)
        .padding(.vertical, 10)
    }
}

private struct AddModulePopover: View {
    @EnvironmentObject private var model: AppModel
    var onPick: (ModuleSpec) -> Void

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            Text("Add a module").font(.headline)
            List {
                ForEach(ModuleGroup.allCases, id: \.self) { grp in
                    Section(grp.title) {
                        ForEach(model.availableModules.filter { $0.group == grp }) { spec in
                            Button {
                                onPick(spec)
                            } label: {
                                HStack {
                                    Text(spec.title)
                                    Spacer()
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
