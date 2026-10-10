// Render anim.html frame by frame in headless Chromium and save PNGs.
//
// usage: node capture.mjs <render dir> <data json> <out dir> [all | f1,f2,...]
import { createRequire } from 'module';
import http from 'http';
import fs from 'fs';
import path from 'path';

const require = createRequire(import.meta.url);
const { chromium } = require('/opt/node22/lib/node_modules/playwright');

const [root, dataFile, outDir, which = 'all'] = process.argv.slice(2);
fs.mkdirSync(outDir, { recursive: true });

const types = { '.js': 'text/javascript', '.mjs': 'text/javascript', '.json': 'application/json',
                '.html': 'text/html', '.woff2': 'font/woff2', '.css': 'text/css' };
const server = http.createServer((req, res) => {
  const file = path.join(root, decodeURIComponent(req.url.split('?')[0]));
  if (!file.startsWith(root)) { res.writeHead(403); res.end(); return; }
  fs.readFile(file, (err, body) => {
    if (err) { res.writeHead(404); res.end(); return; }
    res.writeHead(200, { 'Content-Type': types[path.extname(file)] || 'application/octet-stream' });
    res.end(body);
  });
});
await new Promise((resolve) => server.listen(0, '127.0.0.1', resolve));
const port = server.address().port;

const browser = await chromium.launch({
  args: ['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'],
});
const page = await browser.newPage({ viewport: { width: 1080, height: 1080 }, deviceScaleFactor: 1 });
page.on('pageerror', (e) => console.log('pageerror:', e.message));
page.on('console', (m) => { if (m.type() === 'error') console.log('console:', m.text()); });
await page.goto(`http://127.0.0.1:${port}/${process.env.PAGE || 'anim.html'}?${process.env.QUERY || ('data=' + encodeURIComponent(dataFile))}`);
await page.waitForFunction(() => window.__ready === true, null, { timeout: 120000 });
const total = await page.evaluate(() => window.TOTAL_FRAMES);
const frames = which === 'all' ? [...Array(total).keys()] : which.split(',').map(Number);
console.log(`total frames ${total}; rendering ${frames.length}`);
const t0 = Date.now();
for (const k of frames) {
  await page.evaluate(async (f) => { await window.renderFrame(f); }, k);
  await page.screenshot({ path: path.join(outDir, `${String(k).padStart(4, '0')}.png`) });
  if (k % 30 === 0) console.log(`frame ${k} (${((Date.now() - t0) / 1000).toFixed(1)} s)`);
}
console.log(`done in ${((Date.now() - t0) / 1000).toFixed(1)} s`);
await browser.close();
server.close();
