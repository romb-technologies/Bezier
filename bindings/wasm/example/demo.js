// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
//
// Multi-curve canvas editor for the Bezier WASM bindings. Mirrors the Qt example/:
// create / select / drag / delete curves, with keyboard shortcuts and live visualizations.
// All geometry (sampling, length, extrema, offset, intersections, projection) runs in WASM.
//
// Curve coordinates are in "world" space; a view transform (scale + translate) maps world to
// screen so handles and line widths stay constant size at any zoom.

import { loadBezier, BezierError } from "./dist/index.js";

const { Curve } = await loadBezier();

const stage = document.getElementById("stage");
const canvas = /** @type {HTMLCanvasElement} */ (document.getElementById("canvas"));
const ctx = canvas.getContext("2d");
const readout = document.getElementById("readout");
const help = document.getElementById("help");

const HANDLE_R = 7;
const PICK_DIST = 8;

// --- state ---
let curves = []; // { pts: {x,y}[], selected: boolean }  (world coords)
let show = { box: false, extrema: false, intersections: true, tangent: true };
let offset = 0;
let oHeld = false;
let mode = "idle"; // "placing" | "freehand"
let draft = [];
let cursor = null; // world coords
let drag = null; // { ci, pi }
let freehandDrawing = false;

let view = { scale: 1, tx: 0, ty: 0 };
let panning = false;
let panLast = null;
let DPR = 1;

// --- world <-> screen ---
const w2s = (p) => ({ x: p.x * view.scale + view.tx, y: p.y * view.scale + view.ty });
const s2w = (s) => ({ x: (s.x - view.tx) / view.scale, y: (s.y - view.ty) / view.scale });

// --- helpers ---
const dist2 = (a, b) => (a.x - b.x) ** 2 + (a.y - b.y) ** 2;

function strokePoly(pts) {
  ctx.beginPath();
  pts.forEach((p, k) => {
    const s = w2s(p);
    k ? ctx.lineTo(s.x, s.y) : ctx.moveTo(s.x, s.y);
  });
  ctx.stroke();
}

function dot(p, r, color) {
  const s = w2s(p);
  ctx.beginPath();
  ctx.arc(s.x, s.y, r, 0, Math.PI * 2);
  ctx.fillStyle = color;
  ctx.fill();
}

const selected = () => curves.filter((c) => c.selected);

// --- rendering ---
function draw() {
  const W = canvas.clientWidth;
  const H = canvas.clientHeight;
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0); // device pixels -> CSS px; view applied per-point
  ctx.clearRect(0, 0, W, H);

  const built = curves.map((c) => new Curve(c.pts));
  try {
    curves.forEach((c, k) => {
      const curve = built[k];

      if (show.box) {
        const { min, max } = curve.boundingBox();
        const a = w2s(min);
        const b = w2s(max);
        ctx.strokeStyle = "#e0a800";
        ctx.setLineDash([4, 4]);
        ctx.lineWidth = 1;
        ctx.strokeRect(a.x, a.y, b.x - a.x, b.y - a.y);
        ctx.setLineDash([]);
      }

      if (offset !== 0) {
        const off = Curve.offsetCurve(curve, offset);
        try {
          ctx.strokeStyle = "#2e7d32";
          ctx.setLineDash([6, 4]);
          ctx.lineWidth = 1.5;
          strokePoly(off.polyline());
          ctx.setLineDash([]);
        } finally {
          off.dispose();
        }
      }

      ctx.strokeStyle = c.selected ? "#1565c0" : "#5a5a5a";
      ctx.lineWidth = c.selected ? 2.5 : 2;
      strokePoly(curve.polyline());

      if (show.extrema) {
        for (const t of curve.extrema()) dot(curve.valueAt(t), 4, "#c2185b");
      }

      if (c.selected) {
        ctx.strokeStyle = "#bbb";
        ctx.lineWidth = 1;
        strokePoly(c.pts);
        c.pts.forEach((p) => {
          const s = w2s(p);
          ctx.beginPath();
          ctx.arc(s.x, s.y, HANDLE_R, 0, Math.PI * 2);
          ctx.fillStyle = "#fff";
          ctx.fill();
          ctx.strokeStyle = "#1565c0";
          ctx.lineWidth = 2;
          ctx.stroke();
        });
      }
    });

    if (show.intersections) {
      ctx.fillStyle = "#d32f2f";
      // i === k computes self-intersections of a single curve
      for (let k = 0; k < built.length; k++)
        for (let i = k; i < built.length; i++)
          for (const p of built[k].intersections(built[i])) dot(p, 4, "#d32f2f");
    }

    if (show.tangent && cursor && built.length) {
      let best = null;
      for (const curve of built) {
        const t = curve.projectPoint(cursor);
        const p = curve.valueAt(t);
        const d = dist2(p, cursor);
        if (!best || d < best.d) best = { d, p, tan: curve.tangentAt(t), nor: curve.normalAt(t) };
      }
      const L = 60; // screen px
      const ps = w2s(best.p);
      ctx.strokeStyle = "#aaa";
      ctx.lineWidth = 1;
      strokePoly([cursor, best.p]);
      ctx.strokeStyle = "#1565c0";
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(ps.x - L * best.tan.x, ps.y - L * best.tan.y);
      ctx.lineTo(ps.x + L * best.tan.x, ps.y + L * best.tan.y);
      ctx.stroke();
      ctx.strokeStyle = "#2e7d32";
      ctx.beginPath();
      ctx.moveTo(ps.x, ps.y);
      ctx.lineTo(ps.x + L * best.nor.x, ps.y + L * best.nor.y);
      ctx.stroke();
      dot(best.p, 3, "#1565c0");
    }
  } finally {
    built.forEach((c) => c.dispose());
  }

  if (mode === "placing") {
    if (draft.length >= 2) {
      const preview = new Curve(draft);
      try {
        ctx.strokeStyle = "#1565c0";
        ctx.setLineDash([3, 3]);
        ctx.lineWidth = 2;
        strokePoly(preview.polyline());
        ctx.setLineDash([]);
      } finally {
        preview.dispose();
      }
    }
    draft.forEach((p) => dot(p, 4, "#1565c0"));
  } else if (mode === "freehand" && draft.length >= 2) {
    ctx.strokeStyle = "#888";
    ctx.setLineDash([4, 3]);
    ctx.lineWidth = 1.5;
    strokePoly(draft);
    ctx.setLineDash([]);
  }

  updateReadout();
}

function updateReadout() {
  const zoom = ` · ${Math.round(view.scale * 100)}%`;
  if (mode === "placing") {
    readout.textContent = `placing curve — ${draft.length} point(s); Enter to finish, Esc to cancel`;
    return;
  }
  if (mode === "freehand") {
    readout.textContent = "freehand — drag to draw, release to fit a curve; Esc to cancel";
    return;
  }
  const sel = selected();
  if (sel.length === 1) {
    const c = new Curve(sel[0].pts);
    try {
      readout.textContent = `order ${c.order()} · length ${c.length().toFixed(1)} px${offset ? ` · offset ${offset}` : ""}${zoom}`;
    } finally {
      c.dispose();
    }
  } else {
    readout.textContent = `${curves.length} curve(s), ${sel.length} selected${zoom} · press H for help`;
  }
}

// --- picking (screen coords in) ---
function pointerPos(ev) {
  const r = canvas.getBoundingClientRect();
  return { x: ev.clientX - r.left, y: ev.clientY - r.top };
}

function handleAt(screen) {
  for (let ci = 0; ci < curves.length; ci++) {
    if (!curves[ci].selected) continue;
    const pts = curves[ci].pts;
    for (let pi = 0; pi < pts.length; pi++)
      if (dist2(w2s(pts[pi]), screen) <= (HANDLE_R + 3) ** 2) return { ci, pi };
  }
  return null;
}

function curveAt(screen) {
  const wp = s2w(screen);
  for (let ci = 0; ci < curves.length; ci++) {
    const c = new Curve(curves[ci].pts);
    try {
      const ps = w2s(c.valueAt(c.projectPoint(wp)));
      if (dist2(ps, screen) <= PICK_DIST ** 2) return ci;
    } finally {
      c.dispose();
    }
  }
  return -1;
}

// --- curve ops (round-trip control points through the library) ---
function mutateSelected(fn) {
  for (const c of selected()) {
    const curve = new Curve(c.pts);
    try {
      fn(curve);
      c.pts = curve.controlPoints();
    } catch (e) {
      if (!(e instanceof BezierError)) throw e; // e.g. lowerOrder on a 1st-order curve
    } finally {
      curve.dispose();
    }
  }
  draw();
}

function joinSelected() {
  const sel = selected();
  if (sel.length !== 2) return;
  const a = new Curve(sel[0].pts);
  const b = new Curve(sel[1].pts);
  let pts;
  try {
    const joined = Curve.joinCurves(a, b);
    try {
      pts = joined.controlPoints();
    } finally {
      joined.dispose();
    }
  } finally {
    a.dispose();
    b.dispose();
  }
  curves = curves.filter((c) => !c.selected);
  curves.push({ pts, selected: true });
  draw();
}

function reset() {
  const W = canvas.clientWidth || 800;
  const H = canvas.clientHeight || 500;
  curves = [
    {
      // cubic (order 3)
      pts: [
        { x: W * 0.12, y: H * 0.35 },
        { x: W * 0.35, y: H * 0.78 },
        { x: W * 0.65, y: H * 0.08 },
        { x: W * 0.88, y: H * 0.55 },
      ],
      selected: true,
    },
    {
      // quintic (order 5) — loops over itself and crosses the cubic three times
      pts: [
        { x: W * 0.1, y: H * 0.75 },
        { x: W * 0.55, y: H * 0.05 },
        { x: W * 0.95, y: H * 0.55 },
        { x: W * 0.15, y: H * 0.55 },
        { x: W * 0.55, y: H * 0.95 },
        { x: W * 0.9, y: H * 0.25 },
      ],
      selected: false,
    },
  ];
  mode = "idle";
  draft = [];
  offset = 0;
  view = { scale: 1, tx: 0, ty: 0 };
  draw();
}

// --- mouse ---
canvas.addEventListener("pointerdown", (ev) => {
  const screen = pointerPos(ev);
  const world = s2w(screen);
  if (mode === "placing") {
    draft.push(world);
    draw();
    return;
  }
  if (mode === "freehand") {
    draft = [world];
    freehandDrawing = true;
    canvas.setPointerCapture(ev.pointerId);
    return;
  }
  const h = handleAt(screen);
  if (h) {
    drag = h;
    canvas.setPointerCapture(ev.pointerId);
    return;
  }
  const ci = curveAt(screen);
  if (ci >= 0) {
    if (!ev.shiftKey) curves.forEach((c) => (c.selected = false));
    curves[ci].selected = true;
    draw();
    return;
  }
  // empty space: pan (and clear selection)
  if (!ev.shiftKey) curves.forEach((c) => (c.selected = false));
  panning = true;
  panLast = screen;
  canvas.setPointerCapture(ev.pointerId);
  draw();
});

canvas.addEventListener("pointermove", (ev) => {
  const screen = pointerPos(ev);
  cursor = s2w(screen);
  if (panning) {
    view.tx += screen.x - panLast.x;
    view.ty += screen.y - panLast.y;
    panLast = screen;
    draw();
    return;
  }
  if (freehandDrawing) {
    draft.push(cursor);
    draw();
    return;
  }
  if (drag) {
    curves[drag.ci].pts[drag.pi] = cursor;
    draw();
    return;
  }
  if (show.tangent) draw();
});

canvas.addEventListener("pointerup", () => {
  if (freehandDrawing) {
    freehandDrawing = false;
    if (draft.length >= 2) {
      const fit = Curve.fromPolyline(draft); // auto-order fit
      try {
        curves.forEach((c) => (c.selected = false));
        curves.push({ pts: fit.controlPoints(), selected: true });
      } finally {
        fit.dispose();
      }
    }
    mode = "idle";
    draft = [];
    draw();
    return;
  }
  panning = false;
  drag = null;
});
canvas.addEventListener("pointerleave", () => {
  cursor = null;
  panning = false;
  draw();
});

// zoom centered on the cursor
canvas.addEventListener(
  "wheel",
  (ev) => {
    ev.preventDefault();
    const screen = pointerPos(ev);
    const world = s2w(screen);
    const factor = Math.exp(-ev.deltaY * 0.0015);
    view.scale = Math.min(40, Math.max(0.1, view.scale * factor));
    view.tx = screen.x - world.x * view.scale;
    view.ty = screen.y - world.y * view.scale;
    draw();
  },
  { passive: false },
);

// --- keyboard ---
window.addEventListener("keydown", (ev) => {
  const k = ev.key;
  if (k === "o" || k === "O") oHeld = true;

  if (k === "ArrowUp") {
    ev.preventDefault();
    if (oHeld) {
      offset += 5;
      draw();
    } else mutateSelected((c) => c.raiseOrder());
  } else if (k === "ArrowDown") {
    ev.preventDefault();
    if (oHeld) {
      offset -= 5;
      draw();
    } else mutateSelected((c) => c.lowerOrder());
  } else if (k === "n" || k === "N") {
    mode = "placing";
    draft = [];
    draw();
  } else if (k === "f" || k === "F") {
    mode = "freehand";
    draft = [];
    draw();
  } else if (k === "Enter") {
    if (mode === "placing" && draft.length >= 2) {
      curves.forEach((c) => (c.selected = false));
      curves.push({ pts: draft, selected: true });
      mode = "idle";
      draft = [];
      draw();
    }
  } else if (k === "Escape") {
    mode = "idle";
    draft = [];
    freehandDrawing = false;
    draw();
  } else if (k === "Delete" || k === "Backspace") {
    curves = curves.filter((c) => !c.selected);
    draw();
  } else if (k === "b" || k === "B") {
    show.box = !show.box;
    draw();
  } else if (k === "e" || k === "E") {
    show.extrema = !show.extrema;
    draw();
  } else if (k === "i" || k === "I") {
    show.intersections = !show.intersections;
    draw();
  } else if (k === "p" || k === "P") {
    show.tangent = !show.tangent;
    draw();
  } else if (k === "j" || k === "J") {
    joinSelected();
  } else if (k === "0") {
    view = { scale: 1, tx: 0, ty: 0 };
    draw();
  } else if (k === "h" || k === "H" || k === "?") {
    help.classList.toggle("show");
  } else if (k === "r" || k === "R") {
    reset();
  }
});

window.addEventListener("keyup", (ev) => {
  if (ev.key === "o" || ev.key === "O") oHeld = false;
});

// --- toolbar buttons ---
document.getElementById("new").addEventListener("click", () => {
  mode = "placing";
  draft = [];
  draw();
});
document.getElementById("freehand").addEventListener("click", () => {
  mode = "freehand";
  draft = [];
  draw();
});
document.getElementById("delete").addEventListener("click", () => {
  curves = curves.filter((c) => !c.selected);
  draw();
});
document.getElementById("reset").addEventListener("click", reset);
document.getElementById("help-btn").addEventListener("click", () => help.classList.toggle("show"));

// --- responsive canvas ---
function resize() {
  DPR = window.devicePixelRatio || 1;
  canvas.width = Math.round(stage.clientWidth * DPR);
  canvas.height = Math.round(stage.clientHeight * DPR);
  draw();
}
window.addEventListener("resize", resize);

resize();
reset();
