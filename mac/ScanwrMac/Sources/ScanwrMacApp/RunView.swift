import SwiftUI

struct RunView: View {
    @EnvironmentObject private var model: AppModel

    var body: some View {
        VStack(alignment: .leading, spacing: 12) {
            HStack {
                Text("Run").font(.title2).bold()
                Spacer()
                Button("Run Pipeline") {
                    Task { await model.runPipeline() }
                }
            }

            ScrollView {
                VStack(alignment: .leading, spacing: 6) {
                    ForEach(Array(model.logs.enumerated()), id: \.offset) { _, line in
                        Text(line)
                            .font(.system(.body, design: .monospaced))
                            .frame(maxWidth: .infinity, alignment: .leading)
                    }
                }
            }
            .background(Color(.textBackgroundColor))
            .clipShape(RoundedRectangle(cornerRadius: 8))
        }
        .padding(16)
    }
}

