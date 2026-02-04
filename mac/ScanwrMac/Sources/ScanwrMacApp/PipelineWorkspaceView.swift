import SwiftUI

struct PipelineWorkspaceView: View {
    @EnvironmentObject private var model: AppModel
    @Binding var showSettings: Bool
    @Binding var showConsole: Bool
    @State private var isDropTarget = false

    var body: some View {
        HSplitView {
            PipelineBuilderMainView(isDropTarget: $isDropTarget)

            ModulePaletteSidebar(showSettings: $showSettings, showConsole: $showConsole)
                .frame(minWidth: 280, idealWidth: 320, maxWidth: 360)
        }
    }
}
