import SwiftUI
import UniformTypeIdentifiers

struct PipelineBuilderMainView: View {
    @EnvironmentObject private var model: AppModel
    @Binding var isDropTarget: Bool

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
        .overlay(alignment: .topTrailing) {
            if let selected = model.selectedNodeId,
               let binding = model.nodeBinding(id: selected) {
                NodeInspector(node: binding)
                    .frame(width: 360)
                    .padding(12)
                    .background(.ultraThinMaterial)
                    .clipShape(RoundedRectangle(cornerRadius: 14))
                    .shadow(radius: 10)
                    .padding(24)
            }
        }
    }

    private var header: some View {
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
                guard let id = model.selectedNodeId else { return }
                withAnimation(.snappy(duration: 0.12)) {
                    model.moveNode(id: id, by: -1)
                }
            } label: {
                Label("Up", systemImage: "chevron.up")
            }
            .buttonStyle(.bordered)
            .disabled(model.selectedNodeId == nil || model.nodes.count <= 1)

            Button {
                guard let id = model.selectedNodeId else { return }
                withAnimation(.snappy(duration: 0.12)) {
                    model.moveNode(id: id, by: 1)
                }
            } label: {
                Label("Down", systemImage: "chevron.down")
            }
            .buttonStyle(.bordered)
            .disabled(model.selectedNodeId == nil || model.nodes.count <= 1)

            Button(role: .destructive) {
                model.nodes = []
                model.links = []
                model.selectedNodeId = nil
            } label: {
                Label("Clear", systemImage: "trash")
            }
            .buttonStyle(.bordered)
            .disabled(model.nodes.isEmpty || model.isRunning)
        }
        .padding(12)
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
                List(selection: $model.selectedNodeId) {
                    ForEach(Array(model.nodes.enumerated()), id: \.element.id) { idx, node in
                        PipelineStepRow(stepIndex: idx + 1, node: node)
                            .tag(node.id)
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

    var body: some View {
        let spec = model.spec(for: node.specId)
        let group = spec?.group ?? .pp
        let title = spec?.title ?? node.specId
        let subtitle = spec?.scanpyQualname ?? node.specId

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
        .padding(.vertical, 4)
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
