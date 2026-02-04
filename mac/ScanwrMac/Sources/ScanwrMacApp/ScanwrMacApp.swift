import SwiftUI

@main
struct ScanwrMacApp: App {
    @StateObject private var appModel = AppModel()

    var body: some Scene {
        WindowGroup {
            ShellView()
                .environmentObject(appModel)
        }
        .windowStyle(.titleBar)
    }
}
