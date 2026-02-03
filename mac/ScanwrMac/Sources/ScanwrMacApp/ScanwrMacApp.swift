import SwiftUI

@main
struct ScanwrMacApp: App {
    @StateObject private var appModel = AppModel()

    var body: some Scene {
        WindowGroup {
            ContentView()
                .environmentObject(appModel)
        }
        .windowStyle(.titleBar)
    }
}

