import SwiftUI
import AppKit
import UniformTypeIdentifiers

struct PipelineBuilderMainView: View {
    @EnvironmentObject private var model: AppModel
    @Binding var isDropTarget: Bool
    @State private var draggingNodeId: UUID?

    var body: some View {
        ZStack(alignment: .topLeading) {
            VStack(spacing: 0) {
                header
                Divider()
                content
            }
            .panelChrome()
            .padding(10)

            if isDropTarget {
                Text("Drop to add module")
                    .font(.caption)
                    .padding(.horizontal, 10)
                    .padding(.vertical, 8)
                    .background(.ultraThinMaterial)
                    .clipShape(RoundedRectangle(cornerRadius: 10))
                    .padding(22)
            }
        }
        .onDrop(of: [UTType.plainText], isTargeted: $isDropTarget) { providers in
            handleModuleDrop(providers: providers)
        }
        .onDeleteCommand {
            model.removeSelectedNode()
        }
    }

    private var header: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack(spacing: 10) {
                VStack(alignment: .leading, spacing: 2) {
                    Text("Pipeline")
                        .font(.headline)
                    Text("Double-click a module to append. Drag modules into this panel to add.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }

                Spacer()

                Button {
                    openProjectFolder()
                } label: {
                    Label("Open Folder", systemImage: "folder")
                }
                .buttonStyle(.bordered)
                .disabled(!model.hasProject)

                Button {
                    guard let id = model.selectedNodeId else { return }
                    withAnimation(.snappy(duration: 0.12)) {
                        model.moveNode(id: id, by: -1)
                    }
                } label: {
                    Label("Up", systemImage: "chevron.up")
                }
                .buttonStyle(.bordered)
                .disabled(model.selectedNodeId == nil || model.nodes.count <= 1 || model.isRunning)

                Button {
                    guard let id = model.selectedNodeId else { return }
                    withAnimation(.snappy(duration: 0.12)) {
                        model.moveNode(id: id, by: 1)
                    }
                } label: {
                    Label("Down", systemImage: "chevron.down")
                }
                .buttonStyle(.bordered)
                .disabled(model.selectedNodeId == nil || model.nodes.count <= 1 || model.isRunning)

                Button(role: .destructive) {
                    model.nodes = []
                    model.links = []
                    model.selectedNodeId = nil
                } label: {
                    Label("Clear", systemImage: "trash")
                }
                .buttonStyle(.bordered)
                .disabled(model.nodes.isEmpty || model.isRunning)

                if model.isRunning {
                    Button(role: .destructive) {
                        Task { await model.stopRun() }
                    } label: {
                        Label("Stop", systemImage: "stop.fill")
                    }
                    .buttonStyle(.borderedProminent)
                    .tint(.red)
                    .controlSize(.large)
                } else {
                    Button {
                        model.startRun()
                    } label: {
                        Label("Run Pipeline", systemImage: "play.fill")
                    }
                    .buttonStyle(.borderedProminent)
                    .tint(.green)
                    .controlSize(.large)
                    .disabled(
                        model.outputDirectory.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
                        || model.projectName.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty
                        || model.samples.isEmpty
                        || model.nodes.isEmpty
                    )
                }
            }

            if model.isRunning || model.progressPercent > 0 || !model.progressMessage.isEmpty {
                VStack(alignment: .leading, spacing: 4) {
                    ProgressView(value: model.progressPercent)
                        .progressViewStyle(.linear)
                    if !model.progressMessage.isEmpty {
                        Text(model.progressMessage)
                            .font(.caption)
                            .foregroundStyle(.secondary)
                            .lineLimit(1)
                    }
                }
            }
        }
        .padding(12)
    }

    private func openProjectFolder() {
        guard let url = model.projectPath else { return }
        NSWorkspace.shared.open(url)
    }

    private var content: some View {
        Group {
            if model.nodes.isEmpty {
                VStack(spacing: 10) {
                    Spacer()
                    Text("No modules yet.")
                        .font(.headline)
                    Text("Add modules from the right sidebar to build a step-by-step pipeline.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .multilineTextAlignment(.center)
                        .frame(maxWidth: 420)
                    Spacer()
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .contentShape(Rectangle())
                .onTapGesture { model.selectedNodeId = nil }
            } else {
                List {
                    ForEach(Array(model.nodes.enumerated()), id: \.element.id) { idx, node in
                        PipelineStepRow(
                            stepIndex: idx + 1,
                            node: node,
                            isExpanded: model.selectedNodeId == node.id,
                            draggingNodeId: $draggingNodeId
                        )
                    }
                    .onDelete { indexSet in
                        // Delete supports multi-select; remove by id to avoid index shifting issues.
                        let ids = indexSet.compactMap { i in model.nodes.indices.contains(i) ? model.nodes[i].id : nil }
                        for id in ids { model.removeNode(id: id) }
                    }
                }
                .listStyle(.inset)
            }
        }
        .padding(0)
    }

    private func handleModuleDrop(providers: [NSItemProvider]) -> Bool {
        guard let provider = providers.first else { return false }
        provider.loadItem(forTypeIdentifier: UTType.plainText.identifier, options: nil) { item, _error in
            let text: String? = {
                if let s = item as? String { return s }
                if let data = item as? Data { return String(data: data, encoding: .utf8) }
                if let ns = item as? NSString { return ns as String }
                return nil
            }()
            let specId = (text ?? "").trimmingCharacters(in: .whitespacesAndNewlines)
            guard !specId.isEmpty else { return }

            Task { @MainActor in
                guard let spec = model.availableModules.first(where: { $0.id == specId }) else {
                    model.appendLog("Drop ignored: unknown module id \(specId)")
                    return
                }
                model.appendStep(spec: spec)
            }
        }
        return true
    }
}

private struct PipelineStepRow: View {
    @EnvironmentObject private var model: AppModel
    var stepIndex: Int
    var node: PipelineNode
    var isExpanded: Bool
    @Binding var draggingNodeId: UUID?

    var body: some View {
        let spec = model.spec(for: node.specId)
        let group = spec?.group ?? .pp
        let title = spec?.title ?? node.specId
        let subtitle = spec?.scanpyQualname ?? node.specId

        VStack(alignment: .leading, spacing: 10) {
            HStack(spacing: 12) {
                Text("\(stepIndex)")
                    .font(.caption.weight(.semibold))
                    .foregroundStyle(.secondary)
                    .frame(width: 22, height: 22)
                    .background(Color.primary.opacity(0.08))
                    .clipShape(RoundedRectangle(cornerRadius: 7))

                VStack(alignment: .leading, spacing: 2) {
                    Text(title)
                        .font(.body.weight(.semibold))
                        .lineLimit(1)
                    Text(subtitle)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .lineLimit(1)
                }

                Spacer(minLength: 10)

                HStack(spacing: 6) {
                    Button {
                        withAnimation(.snappy(duration: 0.12)) {
                            model.moveNode(id: node.id, by: -1)
                        }
                    } label: {
                        Image(systemName: "chevron.up")
                    }
                    .buttonStyle(.borderless)
                    .disabled(stepIndex == 1 || model.isRunning)

                    Button {
                        withAnimation(.snappy(duration: 0.12)) {
                            model.moveNode(id: node.id, by: 1)
                        }
                    } label: {
                        Image(systemName: "chevron.down")
                    }
                    .buttonStyle(.borderless)
                    .disabled(stepIndex == model.nodes.count || model.isRunning)
                }
                .foregroundStyle(.secondary)

                Text(group.badge)
                    .font(.caption)
                    .padding(.horizontal, 8)
                    .padding(.vertical, 3)
                    .background(Color(hex: group.colorHex).opacity(0.16))
                    .clipShape(RoundedRectangle(cornerRadius: 8))
            }

            if isExpanded, let binding = model.nodeBinding(id: node.id) {
                NodeInspectorInline(node: binding)
            }
        }
        .padding(.vertical, isExpanded ? 8 : 4)
        .listRowBackground(isExpanded ? Color.primary.opacity(0.05) : Color.clear)
        .contentShape(Rectangle())
        .onTapGesture {
            withAnimation(.snappy(duration: 0.18)) {
                model.selectedNodeId = (model.selectedNodeId == node.id) ? nil : node.id
            }
        }
        .onDrag {
            guard !model.isRunning else { return NSItemProvider() }
            draggingNodeId = node.id
            return NSItemProvider(object: node.id.uuidString as NSString)
        }
        .onDrop(
            of: [UTType.plainText],
            delegate: PipelineStepDropDelegate(
                targetNodeId: node.id,
                model: model,
                draggingNodeId: $draggingNodeId
            )
        )
        .contextMenu {
            Button {
                model.moveNode(id: node.id, by: -1)
            } label: {
                Label("Move Up", systemImage: "chevron.up")
            }
            .disabled(stepIndex == 1)

            Button {
                model.moveNode(id: node.id, by: 1)
            } label: {
                Label("Move Down", systemImage: "chevron.down")
            }
            .disabled(stepIndex == model.nodes.count)

            Button(role: .destructive) {
                model.removeNode(id: node.id)
            } label: {
                Label("Remove", systemImage: "trash")
            }
        }
    }
}

private struct PipelineStepDropDelegate: DropDelegate {
    var targetNodeId: UUID
    var model: AppModel
    @Binding var draggingNodeId: UUID?

    func dropEntered(info: DropInfo) {
        guard !model.isRunning else { return }
        guard let draggingId = draggingNodeId, draggingId != targetNodeId else { return }
        guard let fromIndex = model.nodes.firstIndex(where: { $0.id == draggingId }) else { return }
        guard let toIndex = model.nodes.firstIndex(where: { $0.id == targetNodeId }) else { return }
        guard fromIndex != toIndex else { return }

        withAnimation(.snappy(duration: 0.12)) {
            let destination = (toIndex > fromIndex) ? (toIndex + 1) : toIndex
            model.moveSteps(fromOffsets: IndexSet(integer: fromIndex), toOffset: destination)
        }
    }

    func dropUpdated(info: DropInfo) -> DropProposal? {
        DropProposal(operation: .move)
    }

    func performDrop(info: DropInfo) -> Bool {
        draggingNodeId = nil
        return true
    }
}
