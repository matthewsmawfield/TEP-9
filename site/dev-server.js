#!/usr/bin/env node
// Local preview with live reload, using Node's HTTP server and chokidar.
const chokidar = require('chokidar');
const path = require('path');
const fs = require('fs');
const http = require('http');
const { buildStaticSite } = require('./build.js');
const { HTMLToMarkdownConverter } = require('./html-to-markdown.js');
const TYPES = {'.html':'text/html', '.js':'text/javascript', '.css':'text/css',
    '.json':'application/json', '.png':'image/png', '.jpg':'image/jpeg',
    '.jpeg':'image/jpeg', '.svg':'image/svg+xml', '.pdf':'application/pdf',
    '.woff2':'font/woff2', '.ico':'image/x-icon', '.txt':'text/plain'};

class DevServer {
    constructor() {
        this.isBuilding = false;
        this.buildQueue = false;
        this.clients = new Set();
        this.port = 55536; // Unique port for TEP-9 (Paper 36)
    }
    async build() {
        if (this.isBuilding) { this.buildQueue = true; return; }
        this.isBuilding = true;
        try {
            await buildStaticSite();
            await new HTMLToMarkdownConverter().convertSiteToMarkdown();
            for (const client of this.clients) client.write('data: reload\n\n');
        } finally {
            this.isBuilding = false;
        }
        if (this.buildQueue) { this.buildQueue = false; await this.build(); }
    }
    async start() {
        await this.build();
        const root = fs.realpathSync(path.join(__dirname, 'dist'));
        this.server = http.createServer((req, res) => {
            if (req.method !== 'GET' && req.method !== 'HEAD') {
                res.writeHead(405); return res.end();
            }
            let pathname;
            try { pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname); }
            catch { res.writeHead(400); return res.end(); }
            if (pathname === '/__reload') {
                res.writeHead(200, {'Content-Type':'text/event-stream', 'Cache-Control':'no-cache'});
                this.clients.add(res);
                req.on('close', () => this.clients.delete(res));
                return res.write(': connected\n\n');
            }
            try {
                let file = path.resolve(root, '.' + pathname);
                if (file !== root && !file.startsWith(root + path.sep)) throw Error('outside root');
                if (fs.statSync(file).isDirectory()) file = path.join(file, 'index.html');
                file = fs.realpathSync(file);
                if (!file.startsWith(root + path.sep)) throw Error('outside root');
                const ext = path.extname(file).toLowerCase();
                let content = fs.readFileSync(file);
                if (ext === '.html') content = Buffer.from(content.toString().replace('</body>',
                    '<script>new EventSource("/__reload").onmessage=()=>location.reload()</script></body>'));
                res.writeHead(200, {'Content-Type':TYPES[ext] || 'application/octet-stream',
                    'Cache-Control':'no-store', 'X-Content-Type-Options':'nosniff'});
                res.end(req.method === 'HEAD' ? undefined : content);
            } catch { res.writeHead(404); res.end('Not found'); }
        });
        await new Promise((resolve, reject) => {
            this.server.once('error', reject);
            this.server.listen(this.port, '127.0.0.1', resolve);
        });
        this.watcher = chokidar.watch([
            path.join(__dirname, 'components'), path.join(__dirname, 'index.html'),
            path.join(__dirname, 'manifest.json'), path.join(__dirname, 'claims.json'),
            path.join(__dirname, '..', 'results'),
        ], {ignoreInitial:true});
        this.watcher.on('all', () => this.build().catch(e => console.error(e.message)));
        console.log(`Preview: http://127.0.0.1:${this.port}`);
    }
}
if (require.main === module) new DevServer().start().catch(e => { console.error(e.message); process.exitCode=1; });
module.exports = DevServer;
