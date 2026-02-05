import Foundation
import CoreGraphics

enum ModuleGroup: String, Codable, CaseIterable {
    case pp
    case tl
    case custom

    var title: String {
        switch self {
        case .pp: return "Preprocessing"
        case .tl: return "Tools"
        case .custom: return "Custom"
        }
    }

    var badge: String { rawValue }

    var colorHex: String {
        switch self {
        case .pp: return "#2D6CDF"
        case .tl: return "#1F8A70"
        case .custom: return "#F97316"
        }
    }
}

struct ModuleSpec: Identifiable, Codable, Hashable {
    var id: String
    var group: ModuleGroup
    var namespace: String? // "core" | "experimental" | "external"
    var title: String
    var scanpyQualname: String
}

enum JSONValue: Codable, Hashable {
    case string(String)
    case number(Double)
    case bool(Bool)
    case null

    init(from decoder: Decoder) throws {
        let c = try decoder.singleValueContainer()
        if c.decodeNil() { self = .null; return }
        if let b = try? c.decode(Bool.self) { self = .bool(b); return }
        if let n = try? c.decode(Double.self) { self = .number(n); return }
        if let s = try? c.decode(String.self) { self = .string(s); return }
        throw DecodingError.dataCorruptedError(in: c, debugDescription: "Unsupported JSON value")
    }

    func encode(to encoder: Encoder) throws {
        var c = encoder.singleValueContainer()
        switch self {
        case .string(let s): try c.encode(s)
        case .number(let n): try c.encode(n)
        case .bool(let b): try c.encode(b)
        case .null: try c.encodeNil()
        }
    }

    var stringValue: String? { if case .string(let s) = self { return s }; return nil }
    var boolValue: Bool? { if case .bool(let b) = self { return b }; return nil }
    var doubleValue: Double? { if case .number(let n) = self { return n }; return nil }
}

struct CGPointCodable: Codable, Hashable {
    var x: Double
    var y: Double

    init(_ p: CGPoint) {
        x = Double(p.x)
        y = Double(p.y)
    }

    var cgPoint: CGPoint { CGPoint(x: x, y: y) }
}

struct PipelineNode: Identifiable, Codable, Hashable {
    var id: UUID = UUID()
    var specId: String
    var position: CGPointCodable
    var params: [String: JSONValue]
}

struct PipelineLink: Identifiable, Codable, Hashable {
    var id: UUID = UUID()
    var fromNodeId: UUID
    var toNodeId: UUID
}

struct PipelineRunSummary: Codable, Hashable {
    var outputDir: String
    var results: [SampleRunResult]
}

enum AnalysisMode: String, CaseIterable, Identifiable, Codable, Hashable {
    case concat = "concat"
    case perSample = "per_sample"

    var id: String { rawValue }

    var title: String {
        switch self {
        case .concat: return "Cohort (concat)"
        case .perSample: return "Per-sample"
        }
    }
}

struct SampleMetadata: Identifiable, Codable, Hashable {
    var id: UUID = UUID()
    var sample: String
    var group: String
    var path: String
    // e.g. "scanpy.read_10x_mtx"; empty means "auto"
    var reader: String
}

struct ReaderSuggestion: Codable, Hashable {
    var suggested: String
    var reason: String
}

struct SampleRunResult: Codable, Hashable {
    var sample: String
    var group: String
    var path: String
    var reader: String
    var outputDir: String
    var checkpoints: [String]
    var finalPath: String
    var shape: [Int]
}

struct AdataInspectResult: Codable, Hashable {
    var path: String
    var nObs: Int
    var nVars: Int
    var obsColumns: [String]
    var groupbyCandidates: [String]
    var numericObsColumns: [String]
    var varNames: [String]
    var varNamesTotal: Int
    var varNamesTruncated: Bool
    var layers: [String]
    var hasRaw: Bool
    // Mapping of adata.obsm key -> number of dimensions (columns)
    var obsmDims: [String: Int]
}

struct ViolinPlotRequest: Codable, Hashable {
    var h5adPath: String
    var keys: [String]
    var groupby: String?
    var log: Bool
    var useRaw: Bool?
    var stripplot: Bool
    var jitter: JSONValue?
    var size: Int
    var layer: String?
    var densityNorm: String
    var order: [String]?
    var multiPanel: Bool
    var xlabel: String
    var ylabel: String?
    var rotation: Double?
    var show: Bool?
    var scale: String?
    var kwds: [String: JSONValue]
    var outputPath: String
}

struct ViolinPlotResult: Codable, Hashable {
    var svgPath: String
}

struct CustomPlotRequest: Codable, Hashable {
    // h5ad at .scanwr/checkpoints/<sample>.h5ad
    var h5adPath: String
    // "scatter" | "violin" | "box"
    var plotType: String

    // Key references are encoded strings:
    // - "obs:<column>"
    // - "gene:<var_name>"
    // - "obsm:<key>:<dim>" (dim is 0-based; e.g. obsm:X_umap:0)
    var x: String?
    var y: String?
    var color: String?

    // Expression source
    var layer: String?
    var useRaw: Bool?

    // Generic styling
    var title: String
    var subtitle: String
    var legendTitle: String
    var xLabel: String
    var yLabel: String
    var xTickRotation: Double?

    // Scatter-specific
    var pointSize: Double?
    var alpha: Double?

    // Density-specific
    var densityFill: Bool?

    // Output
    var outputPath: String
}

struct CustomPlotResult: Codable, Hashable {
    var svgPath: String
}
