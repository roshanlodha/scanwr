import SwiftUI

struct RootView: View {
    @EnvironmentObject private var model: AppModel
    @State private var showSettings = false
    @State private var section: SidebarSection = .pipelineBuilder

    var body: some View {
        VStack(spacing: 0) {
            VStack(spacing: 0) {
                NavigationSplitView {
                    Sidebar(section: $section)
                } detail: {
                    switch section {
                    case .metadata:
                        MetadataView()
                    case .pipelineBuilder:
                        PipelineWorkspaceView()
                    case .visualization:
                        VisualizationWorkspaceView(showSettings: $showSettings)
                    case .cohortAnalysis:
                        CohortAnalysisView()
                    }
                }
                .navigationSplitViewStyle(.balanced)

                if model.pipelineBuilderShowConsole && section == .pipelineBuilder {
                    Divider()
                    ConsolePanel(lines: model.logs)
                        .frame(height: 240)
                        .padding(10)
                }
            }
        }
        .task { await model.loadModules() }
        .sheet(isPresented: $showSettings) { SettingsSheet() }
    }
}
