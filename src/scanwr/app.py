from __future__ import annotations

import csv
import os
import tempfile
import tkinter as tk
from dataclasses import dataclass
from tkinter import filedialog, messagebox, ttk
from typing import Dict, List, Optional

from scanwr.pipeline.modules import ModuleInstance, ModuleSpec, available_modules
from scanwr.pipeline.runner import SampleRow, run_pipeline


def _set_safe_env() -> None:
    # Scanpy pulls in matplotlib/fontconfig; make caches land somewhere writable.
    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="scanwr-mpl-"))
    os.environ.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())


@dataclass
class UiSampleRow:
    sample: str
    group: str
    path: str


class ScanwrApp(ttk.Frame):
    def __init__(self, master: tk.Tk):
        super().__init__(master)
        self.master = master
        self.samples: List[UiSampleRow] = []
        self.pipeline: List[ModuleInstance] = []
        self.module_specs: Dict[str, ModuleSpec] = {m.id: m for m in available_modules()}

        self._build_ui()

    def _build_ui(self) -> None:
        self.master.title("scanwr (prototype)")
        self.master.geometry("980x620")

        self.pack(fill="both", expand=True)

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True)

        samples_tab = ttk.Frame(nb)
        pipeline_tab = ttk.Frame(nb)
        nb.add(samples_tab, text="Samples")
        nb.add(pipeline_tab, text="Pipeline Builder")

        self._build_samples_tab(samples_tab)
        self._build_pipeline_tab(pipeline_tab)

    # --------------------
    # Samples tab
    # --------------------
    def _build_samples_tab(self, parent: ttk.Frame) -> None:
        toolbar = ttk.Frame(parent)
        toolbar.pack(fill="x", padx=12, pady=10)

        ttk.Button(toolbar, text="Add Sample", command=self._add_sample_dialog).pack(
            side="left", padx=(0, 8)
        )
        ttk.Button(toolbar, text="Remove Selected", command=self._remove_selected_sample).pack(
            side="left", padx=(0, 8)
        )
        ttk.Button(toolbar, text="Import CSV/TSV", command=self._import_samples).pack(
            side="left", padx=(0, 8)
        )

        cols = ("sample", "group", "path")
        self.samples_tree = ttk.Treeview(parent, columns=cols, show="headings", height=18)
        for c in cols:
            self.samples_tree.heading(c, text=c)
        self.samples_tree.column("sample", width=160, anchor="w")
        self.samples_tree.column("group", width=140, anchor="w")
        self.samples_tree.column("path", width=620, anchor="w")
        self.samples_tree.pack(fill="both", expand=True, padx=12, pady=(0, 12))

    def _refresh_samples_tree(self) -> None:
        for iid in self.samples_tree.get_children():
            self.samples_tree.delete(iid)
        for row in self.samples:
            self.samples_tree.insert("", "end", values=(row.sample, row.group, row.path))

    def _add_sample_dialog(self) -> None:
        win = tk.Toplevel(self.master)
        win.title("Add sample")
        win.resizable(False, False)

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill="both", expand=True)

        sample_var = tk.StringVar()
        group_var = tk.StringVar()
        path_var = tk.StringVar()

        ttk.Label(frm, text="sample").grid(row=0, column=0, sticky="w", pady=(0, 6))
        ttk.Entry(frm, textvariable=sample_var, width=50).grid(row=0, column=1, pady=(0, 6))

        ttk.Label(frm, text="group").grid(row=1, column=0, sticky="w", pady=(0, 6))
        ttk.Entry(frm, textvariable=group_var, width=50).grid(row=1, column=1, pady=(0, 6))

        ttk.Label(frm, text="path").grid(row=2, column=0, sticky="w", pady=(0, 6))
        path_entry = ttk.Entry(frm, textvariable=path_var, width=50)
        path_entry.grid(row=2, column=1, pady=(0, 6))

        def browse() -> None:
            p = filedialog.askdirectory(title="Select input folder") or ""
            if not p:
                p = filedialog.askopenfilename(title="Select input file") or ""
            if p:
                path_var.set(p)

        ttk.Button(frm, text="Browse…", command=browse).grid(row=2, column=2, padx=(8, 0))

        def on_add() -> None:
            sample = sample_var.get().strip()
            group = group_var.get().strip()
            path = path_var.get().strip()
            if not sample or not group or not path:
                messagebox.showerror("Missing fields", "sample, group, and path are required.")
                return
            self.samples.append(UiSampleRow(sample=sample, group=group, path=path))
            self._refresh_samples_tree()
            win.destroy()

        ttk.Button(frm, text="Add", command=on_add).grid(row=3, column=1, sticky="e", pady=(8, 0))
        path_entry.focus_set()

    def _remove_selected_sample(self) -> None:
        sel = self.samples_tree.selection()
        if not sel:
            return
        # Tree is in insertion order; remove by index.
        idxs = sorted((self.samples_tree.index(iid) for iid in sel), reverse=True)
        for idx in idxs:
            if 0 <= idx < len(self.samples):
                self.samples.pop(idx)
        self._refresh_samples_tree()

    def _import_samples(self) -> None:
        p = filedialog.askopenfilename(
            title="Import metadata",
            filetypes=[
                ("CSV", "*.csv"),
                ("TSV", "*.tsv"),
                ("All files", "*.*"),
            ],
        )
        if not p:
            return

        delimiter = "\t" if p.lower().endswith(".tsv") else ","
        try:
            with open(p, newline="") as f:
                reader = csv.DictReader(f, delimiter=delimiter)
                required = {"sample", "group", "path"}
                if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
                    messagebox.showerror(
                        "Invalid file",
                        "Metadata must include columns: sample, group, path",
                    )
                    return
                rows: List[UiSampleRow] = []
                for r in reader:
                    rows.append(
                        UiSampleRow(
                            sample=(r.get("sample") or "").strip(),
                            group=(r.get("group") or "").strip(),
                            path=(r.get("path") or "").strip(),
                        )
                    )
            self.samples = [r for r in rows if r.sample and r.group and r.path]
            self._refresh_samples_tree()
        except Exception as e:
            messagebox.showerror("Import failed", str(e))

    # --------------------
    # Pipeline builder tab
    # --------------------
    def _build_pipeline_tab(self, parent: ttk.Frame) -> None:
        outer = ttk.Frame(parent)
        outer.pack(fill="both", expand=True, padx=12, pady=12)

        outer.columnconfigure(0, weight=1)
        outer.columnconfigure(1, weight=2)
        outer.rowconfigure(0, weight=1)
        outer.rowconfigure(1, weight=0)

        palette_frame = ttk.Labelframe(outer, text="Modules")
        palette_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        pipeline_frame = ttk.Labelframe(outer, text="Pipeline")
        pipeline_frame.grid(row=0, column=1, sticky="nsew")

        self._build_module_palette(palette_frame)
        self._build_pipeline_list(pipeline_frame)

        self.config_frame = ttk.Labelframe(outer, text="Step Settings")
        self.config_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(12, 0))
        self._render_empty_config()

    def _build_module_palette(self, parent: ttk.Frame) -> None:
        self.palette_tree = ttk.Treeview(parent, show="tree", height=18)
        self.palette_tree.pack(fill="both", expand=True, padx=8, pady=8)

        groups = {
            "pp": ("Preprocessing (pp)", "#2D6CDF"),
            "tl": ("Tools (tl)", "#1F8A70"),
            "pl": ("Plotting (pl)", "#B35C00"),
        }
        group_iids: Dict[str, str] = {}
        for g, (label, _) in groups.items():
            group_iids[g] = self.palette_tree.insert("", "end", text=label, open=True)

        for spec in available_modules():
            parent_iid = group_iids.get(spec.group, "")
            self.palette_tree.insert(parent_iid, "end", text=spec.scanpy_qualname, values=(spec.id,))

        add_btn = ttk.Button(parent, text="Add to Pipeline", command=self._add_selected_module)
        add_btn.pack(fill="x", padx=8, pady=(0, 8))

    def _selected_palette_module_id(self) -> Optional[str]:
        sel = self.palette_tree.selection()
        if not sel:
            return None
        iid = sel[0]
        # Leaf nodes have a single value: spec.id
        vals = self.palette_tree.item(iid).get("values") or []
        if not vals:
            return None
        return str(vals[0])

    def _add_selected_module(self) -> None:
        module_id = self._selected_palette_module_id()
        if not module_id:
            return
        spec = self.module_specs[module_id]
        inst = ModuleInstance(spec=spec, params=self._default_params_for(spec))
        self.pipeline.append(inst)
        self._refresh_pipeline_tree()

    def _default_params_for(self, spec: ModuleSpec) -> Dict[str, object]:
        if spec.id == "pp.calculate_qc_metrics":
            return {"qc_vars": "", "percent_top": "50,100", "log1p": True}
        return {}

    def _build_pipeline_list(self, parent: ttk.Frame) -> None:
        cols = ("step", "module")
        self.pipeline_tree = ttk.Treeview(parent, columns=cols, show="headings", height=18)
        self.pipeline_tree.heading("step", text="#")
        self.pipeline_tree.heading("module", text="Module")
        self.pipeline_tree.column("step", width=40, anchor="center")
        self.pipeline_tree.column("module", width=520, anchor="w")
        self.pipeline_tree.pack(fill="both", expand=True, padx=8, pady=8)
        self.pipeline_tree.bind("<<TreeviewSelect>>", lambda _e: self._on_pipeline_select())

        btns = ttk.Frame(parent)
        btns.pack(fill="x", padx=8, pady=(0, 8))
        ttk.Button(btns, text="Remove", command=self._remove_selected_step).pack(
            side="left", padx=(0, 8)
        )
        ttk.Button(btns, text="Move Up", command=lambda: self._move_step(-1)).pack(
            side="left", padx=(0, 8)
        )
        ttk.Button(btns, text="Move Down", command=lambda: self._move_step(1)).pack(
            side="left", padx=(0, 8)
        )

        run_btn = ttk.Button(btns, text="Run", command=self._run_current)
        run_btn.pack(side="right")

        self.log = tk.Text(parent, height=8, wrap="word")
        self.log.pack(fill="both", expand=False, padx=8, pady=(0, 8))
        self.log.configure(state="disabled")

    def _refresh_pipeline_tree(self) -> None:
        for iid in self.pipeline_tree.get_children():
            self.pipeline_tree.delete(iid)
        for i, step in enumerate(self.pipeline, start=1):
            self.pipeline_tree.insert(
                "",
                "end",
                values=(str(i), step.spec.scanpy_qualname),
            )

    def _selected_pipeline_index(self) -> Optional[int]:
        sel = self.pipeline_tree.selection()
        if not sel:
            return None
        return int(self.pipeline_tree.index(sel[0]))

    def _remove_selected_step(self) -> None:
        idx = self._selected_pipeline_index()
        if idx is None:
            return
        if 0 <= idx < len(self.pipeline):
            self.pipeline.pop(idx)
        self._refresh_pipeline_tree()
        self._render_empty_config()

    def _move_step(self, delta: int) -> None:
        idx = self._selected_pipeline_index()
        if idx is None:
            return
        new_idx = idx + delta
        if not (0 <= idx < len(self.pipeline) and 0 <= new_idx < len(self.pipeline)):
            return
        self.pipeline[idx], self.pipeline[new_idx] = self.pipeline[new_idx], self.pipeline[idx]
        self._refresh_pipeline_tree()
        self.pipeline_tree.selection_set(self.pipeline_tree.get_children()[new_idx])

    def _render_empty_config(self) -> None:
        for child in self.config_frame.winfo_children():
            child.destroy()
        ttk.Label(self.config_frame, text="Select a step to configure it.").pack(
            anchor="w", padx=10, pady=10
        )

    def _on_pipeline_select(self) -> None:
        idx = self._selected_pipeline_index()
        if idx is None or not (0 <= idx < len(self.pipeline)):
            self._render_empty_config()
            return
        step = self.pipeline[idx]
        self._render_config_for_step(idx, step)

    def _render_config_for_step(self, idx: int, step: ModuleInstance) -> None:
        for child in self.config_frame.winfo_children():
            child.destroy()

        header = ttk.Frame(self.config_frame)
        header.pack(fill="x", padx=10, pady=(10, 6))
        ttk.Label(header, text=f"Step {idx + 1}: {step.spec.title}").pack(side="left")
        ttk.Label(header, text=step.spec.scanpy_qualname, foreground="#666666").pack(side="right")

        body = ttk.Frame(self.config_frame)
        body.pack(fill="x", padx=10, pady=(0, 10))

        if step.spec.id == "pp.calculate_qc_metrics":
            self._render_calculate_qc_metrics_config(body, step)
            return

        ttk.Label(body, text="No settings available.").pack(anchor="w")

    def _render_calculate_qc_metrics_config(self, parent: ttk.Frame, step: ModuleInstance) -> None:
        qc_vars = tk.StringVar(value=str(step.params.get("qc_vars", "") or ""))
        percent_top = tk.StringVar(value=str(step.params.get("percent_top", "") or ""))
        log1p = tk.BooleanVar(value=bool(step.params.get("log1p", True)))

        def commit() -> None:
            step.params["qc_vars"] = qc_vars.get()
            step.params["percent_top"] = percent_top.get()
            step.params["log1p"] = bool(log1p.get())

        parent.columnconfigure(1, weight=1)

        ttk.Label(parent, text="qc_vars (comma-separated)").grid(row=0, column=0, sticky="w")
        e1 = ttk.Entry(parent, textvariable=qc_vars)
        e1.grid(row=0, column=1, sticky="ew", padx=(8, 0))

        ttk.Label(parent, text="percent_top (e.g. 50,100)").grid(
            row=1, column=0, sticky="w", pady=(8, 0)
        )
        e2 = ttk.Entry(parent, textvariable=percent_top)
        e2.grid(row=1, column=1, sticky="ew", padx=(8, 0), pady=(8, 0))

        cb = ttk.Checkbutton(parent, text="log1p", variable=log1p)
        cb.grid(row=2, column=1, sticky="w", pady=(8, 0))

        ttk.Button(parent, text="Apply", command=commit).grid(row=3, column=1, sticky="e", pady=(10, 0))

    def _append_log(self, text: str) -> None:
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def _run_current(self) -> None:
        if not self.samples:
            messagebox.showerror("No samples", "Add at least one sample first.")
            return
        if not self.pipeline:
            messagebox.showerror("No pipeline", "Add at least one module to the pipeline.")
            return

        samples = [SampleRow(sample=s.sample, group=s.group, path=s.path) for s in self.samples]
        self._append_log(f"Running pipeline with {len(self.pipeline)} step(s) on {len(samples)} sample(s)…")
        try:
            results = run_pipeline(samples=samples, pipeline=self.pipeline)
        except Exception as e:
            self._append_log(f"ERROR: {e}")
            messagebox.showerror("Run failed", str(e))
            return

        for r in results:
            self._append_log(
                f"OK: {r.sample.sample} (group={r.sample.group}) via {r.read.reader_name} "
                f"shape={getattr(r.adata, 'shape', '?')}"
            )
        self._append_log("Done.")


def main() -> int:
    _set_safe_env()

    root = tk.Tk()
    try:
        ttk.Style().theme_use("aqua")
    except tk.TclError:
        pass

    ScanwrApp(root)
    root.mainloop()
    return 0

