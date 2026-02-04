import SwiftUI

struct PipelineWorkspaceView: View {
    @EnvironmentObject private var model: AppModel
    @State private var isDropTarget = false

    var body: some View {
        HSplitView {
            PipelineBuilderMainView(isDropTarget: $isDropTarget)

            ModulePaletteSidebar()
                .frame(minWidth: 280, idealWidth: 320, maxWidth: 360)
        }
    }
}
