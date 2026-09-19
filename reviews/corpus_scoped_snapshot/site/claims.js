// Exact result-leaf bindings for publication numbers. Fail closed on drift.
const fs = require('fs');
const path = require('path');
function renderClaims(html) {
    const registry = JSON.parse(fs.readFileSync(path.join(__dirname, 'claims.json'), 'utf8'));
    return html.replace(/\{\{claim:([a-z0-9_]+)\}\}/g, (_, id) => {
        const spec = registry[id];
        if (!spec) throw new Error(`Unknown claim: ${id}`);
        const file = path.resolve(__dirname, '..', 'results', spec.file);
        if (!file.startsWith(path.resolve(__dirname, '..', 'results') + path.sep)) throw new Error('Invalid claim path');
        let value = JSON.parse(fs.readFileSync(file, 'utf8'));
        for (const key of spec.path) value = value[key];
        if (typeof value !== 'number' || !Number.isFinite(value)) throw new Error(`Invalid numeric leaf: ${id}`);
        if (spec.probability && !(value > 0 && value <= 1)) throw new Error(`Invalid probability: ${id}`);
        const text = spec.format === 'integer' ? String(value) : value.toPrecision(spec.digits || 3);
        return `<span data-claim="${id}">${text}</span>`;
    });
}
module.exports = { renderClaims };
