import * as THREE from "three";
import { OrbitControls } from "https://unpkg.com/three@0.160.0/examples/jsm/controls/OrbitControls.js";
import { STLLoader } from "https://unpkg.com/three@0.160.0/examples/jsm/loaders/STLLoader.js";

const el = Object.fromEntries([
  "surfaceFile", "outputDir", "outputName", "method", "pickerTolerance",
  "resamplingStep", "splineLength", "sectionRadius", "resample", "clipInside",
  "surfaceColor", "centerlineColor", "pointColor", "sectionColor", "surfaceOpacity",
  "showBox", "showAxes", "moveX", "moveY", "moveZ", "rotX", "rotY", "rotZ", "scaleX", "scaleY", "scaleZ",
  "cutMode", "boxXMin", "boxXMax", "boxYMin", "boxYMax", "boxZMin", "boxZMax", "cutResample",
  "pickSource", "pickTarget", "generate", "saveCenterline", "resetPoints",
  "crossSections", "clearSections", "saveCsv", "saveGeometry", "fitBoxToModel", "applyBoxCut", "moveToOrigin", "moveCenterTo", "applyTransform",
  "fitView", "clearLine", "statusText", "viewer", "readout", "log",
].map((id) => [id, document.getElementById(id)]));

const loader = new STLLoader();
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xeef1ef);

const camera = new THREE.PerspectiveCamera(45, 1, 0.001, 100000);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(window.devicePixelRatio || 1);
renderer.localClippingEnabled = false;
el.viewer.prepend(renderer.domElement);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

const raycaster = new THREE.Raycaster();
const pointer = new THREE.Vector2();
const pickMarkers = [null, null];
const selectedPoints = [null, null];

let surfaceMesh = null;
let centerlineObject = null;
let sectionGroup = null;
let boxHelper = null;
let axesHelper = null;
let pickMode = null;
let pointerDown = null;

init();

function init() {
  scene.add(new THREE.AmbientLight(0xffffff, 0.62));
  const key = new THREE.DirectionalLight(0xffffff, 1.35);
  key.position.set(1.8, -2.4, 4.2);
  scene.add(key);
  const fill = new THREE.DirectionalLight(0xffffff, 0.45);
  fill.position.set(-3, 2, -2);
  scene.add(fill);
  axesHelper = new THREE.AxesHelper(25);
  scene.add(axesHelper);

  el.outputDir.value = "";
  el.surfaceFile.addEventListener("change", loadSelectedFile);
  el.pickSource.addEventListener("click", () => setPickMode("source"));
  el.pickTarget.addEventListener("click", () => setPickMode("target"));
  el.generate.addEventListener("click", generateCenterline);
  el.saveCenterline.addEventListener("click", saveCenterline);
  el.crossSections.addEventListener("click", calculateCrossSections);
  el.clearSections.addEventListener("click", clearCrossSections);
  el.saveCsv.addEventListener("click", saveCsv);
  el.saveGeometry.addEventListener("click", saveGeometry);
  el.resetPoints.addEventListener("click", resetPoints);
  el.fitBoxToModel.addEventListener("click", fitBoxToModel);
  el.applyBoxCut.addEventListener("click", applyBoxCut);
  el.moveToOrigin.addEventListener("click", () => applyTransform("origin"));
  el.moveCenterTo.addEventListener("click", () => applyTransform("center"));
  el.applyTransform.addEventListener("click", () => applyTransform("apply"));
  el.fitView.addEventListener("click", fitView);
  el.clearLine.addEventListener("click", clearCenterline);
  el.showBox.addEventListener("change", updateBoxVisibility);
  el.showAxes.addEventListener("change", updateAxesVisibility);
  el.surfaceColor.addEventListener("input", applyDisplay);
  el.centerlineColor.addEventListener("input", applyDisplay);
  el.pointColor.addEventListener("input", applyDisplay);
  el.sectionColor.addEventListener("input", applyDisplay);
  el.surfaceOpacity.addEventListener("input", applyDisplay);
  for (const input of boxInputs()) input.addEventListener("input", updateBoxPreview);
  renderer.domElement.addEventListener("pointerdown", onPointerDown);
  renderer.domElement.addEventListener("pointerup", onPointerUp);
  window.addEventListener("resize", resize);
  resize();
  animate();
  log("Ready. Load an STL file to begin.");
  updateAxesVisibility();
}

async function loadSelectedFile() {
  const file = el.surfaceFile.files[0];
  if (!file) return;
  clearScene();
  el.statusText.textContent = "Loading STL...";
  try {
    const buffer = await file.arrayBuffer();
    const geometry = loader.parse(buffer);
    geometry.computeVertexNormals();
    geometry.computeBoundingBox();

    surfaceMesh = new THREE.Mesh(
      geometry,
      new THREE.MeshStandardMaterial({
        color: el.surfaceColor.value,
        roughness: 0.72,
        metalness: 0.02,
        transparent: true,
        opacity: Number(el.surfaceOpacity.value),
        side: THREE.DoubleSide,
      }),
    );
    scene.add(surfaceMesh);

    await uploadModel(file);
    el.outputName.value = file.name.replace(/\.[^.]+$/, "") || "centerline";
    updateAxesScale();
    fitBoxToModel();
    fitView();
    el.statusText.textContent = `${file.name} loaded`;
    log(`Loaded ${file.name}`);
  } catch (error) {
    el.statusText.textContent = "Load failed";
    logError(error);
  }
}

async function uploadModel(file) {
  const res = await fetch("/api/model", {
    method: "POST",
    headers: { "X-Filename": encodeURIComponent(file.name) },
    body: file,
  });
  const payload = await readJson(res);
  log(`Backend model: ${payload.points} points, ${payload.cells} cells`);
}

function setPickMode(mode) {
  if (!surfaceMesh) {
    log("Load an STL before picking points.");
    return;
  }
  pickMode = pickMode === mode ? null : mode;
  el.pickSource.classList.toggle("active", pickMode === "source");
  el.pickTarget.classList.toggle("active", pickMode === "target");
  el.statusText.textContent = pickMode ? `Click the surface to pick ${pickMode}` : "Camera interaction";
}

function onPointerDown(event) {
  pointerDown = { x: event.clientX, y: event.clientY };
}

function onPointerUp(event) {
  if (!pickMode || !surfaceMesh || !pointerDown) return;
  const dx = event.clientX - pointerDown.x;
  const dy = event.clientY - pointerDown.y;
  if (Math.hypot(dx, dy) > 4) return;

  const hit = pickSurface(event);
  if (!hit) {
    log("No valid surface point selected.");
    return;
  }
  const index = pickMode === "source" ? 0 : 1;
  selectedPoints[index] = [hit.point.x, hit.point.y, hit.point.z];
  setMarker(index, hit.point);
  log(`Selected ${pickMode}: ${formatPoint(selectedPoints[index])}`);
  pickMode = null;
  el.pickSource.classList.remove("active");
  el.pickTarget.classList.remove("active");
  el.statusText.textContent = "Camera interaction";
  updateReadout();
}

function pickSurface(event) {
  const rect = renderer.domElement.getBoundingClientRect();
  pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
  raycaster.setFromCamera(pointer, camera);
  const hits = raycaster.intersectObject(surfaceMesh, false);
  return hits[0] || null;
}

function setMarker(index, position) {
  if (pickMarkers[index]) scene.remove(pickMarkers[index]);
  const radius = markerRadius();
  const marker = new THREE.Mesh(
    new THREE.SphereGeometry(radius, 24, 12),
    new THREE.MeshStandardMaterial({ color: el.pointColor.value, roughness: 0.5 }),
  );
  marker.position.copy(position);
  pickMarkers[index] = marker;
  scene.add(marker);
}

async function generateCenterline() {
  if (!selectedPoints[0] || !selectedPoints[1]) {
    log("Pick source and target before generating centerline.");
    return;
  }
  el.statusText.textContent = "Generating centerline...";
  try {
    const payload = {
      source: selectedPoints[0],
      target: selectedPoints[1],
      outputDir: el.outputDir.value,
      outputName: el.outputName.value,
      config: readConfig(),
    };
    log(`Generating: ${JSON.stringify(payload.config)}`);
    const res = await fetch("/api/centerline", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await readJson(res);
    drawCenterline(data.points);
    el.statusText.textContent = `Centerline generated (${data.count} points)`;
    log(`Centerline generated: ${data.count} points`);
  } catch (error) {
    el.statusText.textContent = "Generation failed";
    logError(error);
  }
}

async function saveCenterline() {
  try {
    const res = await fetch("/api/save-centerline", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ outputDir: el.outputDir.value, outputName: el.outputName.value }),
    });
    const data = await readJson(res);
    log(`Centerline saved: ${data.path}`);
  } catch (error) {
    logError(error);
  }
}

function drawCenterline(points) {
  clearCenterline();
  const geometry = new THREE.BufferGeometry().setFromPoints(points.map((p) => new THREE.Vector3(p[0], p[1], p[2])));
  const line = new THREE.Line(
    geometry,
    new THREE.LineBasicMaterial({ color: el.centerlineColor.value }),
  );
  centerlineObject = line;
  scene.add(centerlineObject);
}

async function calculateCrossSections() {
  el.statusText.textContent = "Calculating cross sections...";
  try {
    const res = await fetch("/api/cross-sections", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        outputDir: el.outputDir.value,
        outputName: el.outputName.value,
        config: readConfig(),
      }),
    });
    const data = await readJson(res);
    drawCrossSections(data.sections);
    el.statusText.textContent = `Cross sections calculated (${data.count})`;
    log(`Cross sections: ${data.count}; first area=${data.areas[0]?.toFixed?.(6) ?? "n/a"}`);
  } catch (error) {
    el.statusText.textContent = "Cross section failed";
    logError(error);
  }
}

async function saveCsv() {
  try {
    const res = await fetch("/api/save-csv", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ outputDir: el.outputDir.value, outputName: el.outputName.value }),
    });
    const data = await readJson(res);
    log(`CSV saved: ${data.path}`);
  } catch (error) {
    logError(error);
  }
}

function drawCrossSections(sections) {
  clearCrossSections();
  sectionGroup = new THREE.Group();
  const material = new THREE.LineBasicMaterial({ color: el.sectionColor.value, transparent: true, opacity: 0.78 });
  for (const loops of sections) {
    for (const points of loops) {
      if (!points || points.length < 2) continue;
      const vertices = points.map((p) => new THREE.Vector3(p[0], p[1], p[2]));
      if (!vertices[0].equals(vertices[vertices.length - 1])) vertices.push(vertices[0].clone());
      const geometry = new THREE.BufferGeometry().setFromPoints(vertices);
      sectionGroup.add(new THREE.Line(geometry, material));
    }
  }
  scene.add(sectionGroup);
}

function clearCrossSections() {
  if (sectionGroup) scene.remove(sectionGroup);
  sectionGroup = null;
}

function resetPoints() {
  selectedPoints[0] = null;
  selectedPoints[1] = null;
  for (let i = 0; i < pickMarkers.length; i += 1) {
    if (pickMarkers[i]) scene.remove(pickMarkers[i]);
    pickMarkers[i] = null;
  }
  updateReadout();
  log("Selection reset.");
}

function clearCenterline() {
  if (centerlineObject) scene.remove(centerlineObject);
  centerlineObject = null;
  clearCrossSections();
}

function clearScene() {
  if (surfaceMesh) scene.remove(surfaceMesh);
  surfaceMesh = null;
  clearCenterline();
  resetPoints();
  if (boxHelper) scene.remove(boxHelper);
  boxHelper = null;
}

function readConfig() {
  return {
    method: el.method.value,
    pickerTolerance: number(el.pickerTolerance, 0.005),
    resamplingStep: number(el.resamplingStep, 0.05),
    splineLength: number(el.splineLength, 0.5),
    sectionRadius: number(el.sectionRadius, 20),
    resample: el.resample.checked,
    clipInside: el.clipInside.checked,
  };
}

function applyDisplay() {
  if (surfaceMesh) {
    surfaceMesh.material.color.set(el.surfaceColor.value);
    surfaceMesh.material.opacity = Number(el.surfaceOpacity.value);
    surfaceMesh.material.needsUpdate = true;
  }
  if (centerlineObject) centerlineObject.material.color.set(el.centerlineColor.value);
  if (sectionGroup) {
    sectionGroup.children.forEach((child) => child.material.color.set(el.sectionColor.value));
  }
  for (const marker of pickMarkers) {
    if (marker) marker.material.color.set(el.pointColor.value);
  }
}

function fitBoxToModel() {
  if (!surfaceMesh) return;
  surfaceMesh.geometry.computeBoundingBox();
  const b = surfaceMesh.geometry.boundingBox;
  const pad = Math.max(b.max.x - b.min.x, b.max.y - b.min.y, b.max.z - b.min.z) * 0.04;
  el.boxXMin.value = cleanNumber(b.min.x + pad);
  el.boxXMax.value = cleanNumber(b.max.x - pad);
  el.boxYMin.value = cleanNumber(b.min.y + pad);
  el.boxYMax.value = cleanNumber(b.max.y - pad);
  el.boxZMin.value = cleanNumber(b.min.z + pad);
  el.boxZMax.value = cleanNumber(b.max.z - pad);
  updateBoxPreview();
}

function updateBoxPreview() {
  const bounds = readBoxBounds(false);
  if (!bounds) return;
  const box = new THREE.Box3(
    new THREE.Vector3(bounds[0], bounds[2], bounds[4]),
    new THREE.Vector3(bounds[1], bounds[3], bounds[5]),
  );
  if (boxHelper) scene.remove(boxHelper);
  boxHelper = makeBoxObject(box);
  boxHelper.visible = el.showBox.checked;
  scene.add(boxHelper);
}

function updateBoxVisibility() {
  if (boxHelper) boxHelper.visible = el.showBox.checked;
}

function updateAxesVisibility() {
  if (axesHelper) axesHelper.visible = el.showAxes.checked;
}

async function applyBoxCut() {
  if (!surfaceMesh) {
    log("Load an STL before cutting geometry.");
    return;
  }
  const bounds = readBoxBounds(true);
  if (!bounds) return;
  el.statusText.textContent = "Cutting geometry...";
  try {
    const res = await fetch("/api/cut-box", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        bounds,
        mode: el.cutMode.value === "inside" ? "inside" : "outside",
        resample: el.cutResample.checked,
      }),
    });
    const data = await readJson(res);
    await reloadCurrentSurface();
    clearCenterline();
    resetPoints();
    fitBoxToModel();
    fitView();
    el.statusText.textContent = `Geometry cut (${data.points} points)`;
    log(`Geometry cut: ${data.points} points, ${data.cells} cells`);
  } catch (error) {
    el.statusText.textContent = "Cut failed";
    logError(error);
  }
}

async function applyTransform(action) {
  if (!surfaceMesh) {
    log("Load an STL before transforming geometry.");
    return;
  }
  el.statusText.textContent = "Transforming geometry...";
  try {
    const payload = {
      action,
      translate: [number(el.moveX, 0), number(el.moveY, 0), number(el.moveZ, 0)],
      rotate: [number(el.rotX, 0), number(el.rotY, 0), number(el.rotZ, 0)],
      scale: [number(el.scaleX, 1), number(el.scaleY, 1), number(el.scaleZ, 1)],
    };
    const res = await fetch("/api/transform", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await readJson(res);
    await reloadCurrentSurface();
    clearCenterline();
    resetPoints();
    fitBoxToModel();
    fitView();
    el.statusText.textContent = `Geometry transformed (${data.points} points)`;
    log(`Geometry transformed: ${data.points} points, bounds=${data.bounds.map((v) => Number(v).toFixed(3)).join(", ")}`);
  } catch (error) {
    el.statusText.textContent = "Transform failed";
    logError(error);
  }
}

async function saveGeometry() {
  try {
    const res = await fetch("/api/save-stl", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ outputDir: el.outputDir.value, outputName: el.outputName.value }),
    });
    const data = await readJson(res);
    log(`Edited STL saved: ${data.path}`);
  } catch (error) {
    logError(error);
  }
}

async function reloadCurrentSurface() {
  const res = await fetch(`/api/surface.stl?t=${Date.now()}`);
  if (!res.ok) throw new Error(res.statusText);
  const buffer = await res.arrayBuffer();
  const geometry = loader.parse(buffer);
  geometry.computeVertexNormals();
  geometry.computeBoundingBox();
  if (surfaceMesh) scene.remove(surfaceMesh);
  surfaceMesh = new THREE.Mesh(
    geometry,
    new THREE.MeshStandardMaterial({
      color: el.surfaceColor.value,
      roughness: 0.72,
      metalness: 0.02,
      transparent: true,
      opacity: Number(el.surfaceOpacity.value),
      side: THREE.DoubleSide,
    }),
  );
  scene.add(surfaceMesh);
  updateAxesScale();
}

function makeBoxObject(box) {
  const group = new THREE.Group();
  const size = box.getSize(new THREE.Vector3());
  const center = box.getCenter(new THREE.Vector3());
  const fill = new THREE.Mesh(
    new THREE.BoxGeometry(size.x, size.y, size.z),
    new THREE.MeshBasicMaterial({ color: 0x2f7d68, transparent: true, opacity: 0.13, depthWrite: false }),
  );
  fill.position.copy(center);
  const edges = new THREE.LineSegments(
    new THREE.EdgesGeometry(fill.geometry),
    new THREE.LineBasicMaterial({ color: 0x1a786f, transparent: true, opacity: 0.95 }),
  );
  edges.position.copy(center);
  group.add(fill, edges);
  return group;
}

function readBoxBounds(showErrors) {
  const bounds = [
    number(el.boxXMin, NaN), number(el.boxXMax, NaN),
    number(el.boxYMin, NaN), number(el.boxYMax, NaN),
    number(el.boxZMin, NaN), number(el.boxZMax, NaN),
  ];
  const invalid = bounds.some((v) => !Number.isFinite(v)) || bounds[0] >= bounds[1] || bounds[2] >= bounds[3] || bounds[4] >= bounds[5];
  if (invalid) {
    if (showErrors) log("Invalid box range: min values must be smaller than max values.");
    return null;
  }
  return bounds;
}

function boxInputs() {
  return [el.boxXMin, el.boxXMax, el.boxYMin, el.boxYMax, el.boxZMin, el.boxZMax];
}

function fitView() {
  if (!surfaceMesh) return;
  surfaceMesh.geometry.computeBoundingBox();
  const box = surfaceMesh.geometry.boundingBox.clone();
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3());
  const radius = Math.max(size.x, size.y, size.z) || 1;
  controls.target.copy(center);
  camera.near = Math.max(radius / 10000, 0.001);
  camera.far = radius * 10000;
  camera.position.copy(center).add(new THREE.Vector3(radius * 1.1, -radius * 1.8, radius * 0.9));
  camera.updateProjectionMatrix();
  controls.update();
}

function updateAxesScale() {
  if (!surfaceMesh || !axesHelper) return;
  surfaceMesh.geometry.computeBoundingBox();
  const box = surfaceMesh.geometry.boundingBox;
  const size = box.getSize(new THREE.Vector3());
  const radius = Math.max(size.x, size.y, size.z) || 1;
  scene.remove(axesHelper);
  axesHelper = new THREE.AxesHelper(radius * 0.25);
  axesHelper.visible = el.showAxes.checked;
  scene.add(axesHelper);
}

function markerRadius() {
  if (!surfaceMesh) return 1;
  const box = surfaceMesh.geometry.boundingBox || new THREE.Box3().setFromObject(surfaceMesh);
  const size = box.getSize(new THREE.Vector3());
  return Math.max(size.x, size.y, size.z) * 0.006 || 1;
}

function resize() {
  const rect = el.viewer.getBoundingClientRect();
  const width = Math.max(1, Math.floor(rect.width));
  const height = Math.max(1, Math.floor(rect.height));
  renderer.setSize(width, height, false);
  camera.aspect = width / height;
  camera.updateProjectionMatrix();
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}

async function readJson(res) {
  const data = await res.json();
  if (!res.ok || data.error) throw new Error(data.error || res.statusText);
  return data;
}

function number(input, fallback) {
  const value = Number(input.value);
  return Number.isFinite(value) ? value : fallback;
}

function formatPoint(point) {
  return point.map((v) => v.toFixed(6)).join(", ");
}

function cleanNumber(value) {
  return Number(value).toFixed(6).replace(/\.?0+$/, "");
}

function updateReadout() {
  el.readout.innerHTML = `source: ${selectedPoints[0] ? formatPoint(selectedPoints[0]) : "none"}<br>target: ${selectedPoints[1] ? formatPoint(selectedPoints[1]) : "none"}`;
}

function log(message) {
  el.log.textContent += `${message}\n`;
  el.log.scrollTop = el.log.scrollHeight;
}

function logError(error) {
  log(`Error: ${error.message || error}`);
}
