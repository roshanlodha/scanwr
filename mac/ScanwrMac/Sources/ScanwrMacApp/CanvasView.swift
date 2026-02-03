import SwiftUI

struct CanvasView: View {
    @EnvironmentObject private var model: AppModel

    @State private var linkingFromNodeId: UUID?
    @State private var linkingDragPoint: CGPoint = .zero

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .topLeading) {
                DotGridBackground()

                LinksLayer(
                    links: model.links,
                    nodes: model.nodes,
                    linkingFrom: linkingFromNodeId,
                    linkingDragPoint: linkingDragPoint
                )

                ForEach(model.nodes) { node in
                    NodeBubble(
                        node: node,
                        isSelected: model.selectedNodeId == node.id,
                        spec: model.spec(for: node.specId),
                        onSelect: {
                            model.selectedNodeId = node.id
                        },
                        onDrag: { delta in
                            if let idx = model.nodes.firstIndex(where: { $0.id == node.id }) {
                                var n = model.nodes[idx]
                                let p = n.position.cgPoint
                                n.position = CGPointCodable(CGPoint(x: p.x + delta.width, y: p.y + delta.height))
                                model.nodes[idx] = n
                            }
                        },
                        onStartLink: { startPoint in
                            linkingFromNodeId = node.id
                            linkingDragPoint = startPoint
                        },
                        onUpdateLink: { p in
                            linkingDragPoint = p
                        },
                        onEndLink: { dropPoint in
                            defer { linkingFromNodeId = nil }
                            guard let from = linkingFromNodeId else { return }
                            if let target = nearestInputNode(at: dropPoint, excluding: from) {
                                model.addLink(from: from, to: target)
                            }
                        }
                    )
                    .position(node.position.cgPoint)
                }

                if let selected = model.selectedNodeId,
                   let binding = model.nodeBinding(id: selected) {
                    NodeInspector(node: binding)
                        .frame(width: 360)
                        .padding(12)
                        .background(.ultraThinMaterial)
                        .clipShape(RoundedRectangle(cornerRadius: 14))
                        .shadow(radius: 10)
                        .padding(14)
                        .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .topTrailing)
                }
            }
            .contentShape(Rectangle())
            .onTapGesture {
                model.selectedNodeId = nil
            }
            .frame(width: geo.size.width, height: geo.size.height)
        }
    }

    private func nearestInputNode(at point: CGPoint, excluding: UUID) -> UUID? {
        // Simple hit-test: pick the nearest node center within a radius.
        let radius: CGFloat = 70
        var best: (id: UUID, d: CGFloat)?
        for n in model.nodes where n.id != excluding {
            let c = n.position.cgPoint
            let dx = c.x - point.x
            let dy = c.y - point.y
            let d = sqrt(dx * dx + dy * dy)
            if d <= radius && (best == nil || d < best!.d) {
                best = (n.id, d)
            }
        }
        return best?.id
    }
}

private struct DotGridBackground: View {
    var body: some View {
        Canvas { ctx, size in
            let spacing: CGFloat = 18
            let dotRadius: CGFloat = 1.2
            let color = Color.primary.opacity(0.08)

            var path = Path()
            var y: CGFloat = 0
            while y <= size.height {
                var x: CGFloat = 0
                while x <= size.width {
                    path.addEllipse(in: CGRect(x: x - dotRadius, y: y - dotRadius, width: dotRadius * 2, height: dotRadius * 2))
                    x += spacing
                }
                y += spacing
            }
            ctx.fill(path, with: .color(color))
        }
        .background(Color(NSColor.windowBackgroundColor))
    }
}

private struct LinksLayer: View {
    var links: [PipelineLink]
    var nodes: [PipelineNode]
    var linkingFrom: UUID?
    var linkingDragPoint: CGPoint

    var body: some View {
        Canvas { ctx, _ in
            func center(_ id: UUID) -> CGPoint? {
                nodes.first(where: { $0.id == id })?.position.cgPoint
            }

            for l in links {
                guard let a = center(l.fromNodeId), let b = center(l.toNodeId) else { continue }
                let p = curvedPath(from: a, to: b)
                ctx.stroke(p, with: .color(Color.primary.opacity(0.35)), lineWidth: 3)
            }

            if let from = linkingFrom, let a = center(from) as CGPoint? {
                let p = curvedPath(from: a, to: linkingDragPoint)
                ctx.stroke(p, with: .color(Color.accentColor.opacity(0.7)), style: StrokeStyle(lineWidth: 3, dash: [8, 6]))
            }
        }
    }

    private func curvedPath(from: CGPoint, to: CGPoint) -> Path {
        var path = Path()
        path.move(to: from)
        let dx = max(40, abs(to.x - from.x) * 0.45)
        let c1 = CGPoint(x: from.x + dx, y: from.y)
        let c2 = CGPoint(x: to.x - dx, y: to.y)
        path.addCurve(to: to, control1: c1, control2: c2)
        return path
    }
}

private struct NodeBubble: View {
    var node: PipelineNode
    var isSelected: Bool
    var spec: ModuleSpec?

    var onSelect: () -> Void
    var onDrag: (CGSize) -> Void
    var onStartLink: (CGPoint) -> Void
    var onUpdateLink: (CGPoint) -> Void
    var onEndLink: (CGPoint) -> Void

    @State private var dragAccum: CGSize = .zero

    var body: some View {
        let group = spec?.group ?? .pp
        let title = spec?.title ?? node.specId
        let badge = group.badge

        ZStack {
            RoundedRectangle(cornerRadius: 999)
                .fill(Color(hex: group.colorHex).opacity(isSelected ? 0.24 : 0.16))
                .overlay(
                    RoundedRectangle(cornerRadius: 999)
                        .stroke(isSelected ? Color(hex: group.colorHex).opacity(0.9) : Color.primary.opacity(0.15), lineWidth: 2)
                )

            HStack(spacing: 10) {
                Circle()
                    .fill(Color.primary.opacity(0.18))
                    .frame(width: 12, height: 12)
                    .overlay(Circle().stroke(Color.primary.opacity(0.18), lineWidth: 1))

                VStack(alignment: .leading, spacing: 2) {
                    Text(title)
                        .font(.headline)
                        .lineLimit(1)
                    Text(spec?.scanpyQualname ?? node.specId)
                        .font(.caption)
                        .foregroundStyle(.secondary)
                        .lineLimit(1)
                }
                Spacer(minLength: 10)
                Text(badge)
                    .font(.caption)
                    .padding(.horizontal, 8)
                    .padding(.vertical, 4)
                    .background(Color(hex: group.colorHex).opacity(0.18))
                    .clipShape(RoundedRectangle(cornerRadius: 10))

                OutputPort(onStart: onStartLink, onUpdate: onUpdateLink, onEnd: onEndLink)
            }
            .padding(.horizontal, 14)
            .padding(.vertical, 10)
        }
        .frame(width: 420, height: 78)
        .onTapGesture {
            onSelect()
        }
        .gesture(
            DragGesture(minimumDistance: 2)
                .onChanged { v in
                    let delta = CGSize(width: v.translation.width - dragAccum.width, height: v.translation.height - dragAccum.height)
                    dragAccum = v.translation
                    onDrag(delta)
                }
                .onEnded { _ in
                    dragAccum = .zero
                }
        )
    }
}

private struct OutputPort: View {
    var onStart: (CGPoint) -> Void
    var onUpdate: (CGPoint) -> Void
    var onEnd: (CGPoint) -> Void

    var body: some View {
        GeometryReader { geo in
            Circle()
                .fill(Color.accentColor.opacity(0.9))
                .overlay(Circle().stroke(Color.white.opacity(0.5), lineWidth: 1))
                .frame(width: 14, height: 14)
                .contentShape(Circle())
                .gesture(
                    DragGesture(minimumDistance: 0)
                        .onChanged { v in
                            let p = CGPoint(x: geo.frame(in: .global).midX + v.translation.width,
                                            y: geo.frame(in: .global).midY + v.translation.height)
                            if v.startLocation == v.location {
                                onStart(p)
                            } else {
                                onUpdate(p)
                            }
                        }
                        .onEnded { v in
                            let p = CGPoint(x: geo.frame(in: .global).midX + v.translation.width,
                                            y: geo.frame(in: .global).midY + v.translation.height)
                            onEnd(p)
                        }
                )
        }
        .frame(width: 16, height: 16)
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

