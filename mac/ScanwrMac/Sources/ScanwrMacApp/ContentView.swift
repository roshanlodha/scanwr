import SwiftUI

struct ContentView: View {
    @EnvironmentObject private var model: AppModel

    var body: some View {
        NavigationSplitView {
            List {
                NavigationLink("Samples") { SamplesView() }
                NavigationLink("Pipeline Builder") { PipelineBuilderView() }
                NavigationLink("Run") { RunView() }
            }
            .frame(minWidth: 180)
        } detail: {
            SamplesView()
        }
        .task {
            await model.loadModules()
        }
    }
}

