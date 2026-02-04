import SwiftUI

struct ShellView: View {
    @EnvironmentObject private var model: AppModel

    var body: some View {
        if model.hasProject {
            RootView()
        } else {
            WelcomeView()
        }
    }
}

