// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "ScanwrMac",
    platforms: [
        // Keep reasonably recent; adjust upward once you confirm the exact macOS version naming you want.
        .macOS(.v14),
    ],
    products: [
        .executable(name: "ScanwrMacApp", targets: ["ScanwrMacApp"])
    ],
    dependencies: [],
    targets: [
        .executableTarget(
            name: "ScanwrMacApp",
            resources: [
                .process("Resources")
            ],
            swiftSettings: [
                // Prototype: relax strict concurrency so we can iterate quickly.
                // We can tighten this back up once the architecture stabilizes.
                .unsafeFlags(["-Xfrontend", "-strict-concurrency=minimal"]),
            ]
        )
    ],
    swiftLanguageModes: [.v5]
)
