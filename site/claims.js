// Exact result-leaf bindings for publication numbers. Fail closed on drift.
const fs = require('fs');
const path = require('path');
function renderClaims(html) {
    const registry = JSON.parse(fs.readFileSync(path.join(__dirname, 'claims.json'), 'utf8'));
    return html.replace(/\{\{claim:([a-z0-9_]+)\}\}/g, (match, id, offset) => {
        const spec = registry[id];
        if (!spec) throw new Error(`Unknown claim: ${id}`);
        const file = path.resolve(__dirname, '..', 'results', spec.file);
        if (!file.startsWith(path.resolve(__dirname, '..', 'results') + path.sep)) throw new Error('Invalid claim path');
        let value = JSON.parse(fs.readFileSync(file, 'utf8'));
        for (const key of spec.path) value = value[key];
        if (typeof value !== 'number' || !Number.isFinite(value)) throw new Error(`Invalid numeric leaf: ${id}`);
        if (spec.probability && !(value > 0 && value <= 1)) throw new Error(`Invalid probability: ${id}`);
        // A claim inside $...$ must emit plain text (or TeX for sci format):
        // an HTML element would split the formula and corrupt MathJax parsing.
        const inMath = html.slice(0, offset).replace(/<[a-zA-Z/!][^>]*>/g, '').split('$').length % 2 === 0;
        let text;
        if (inMath && spec.format === 'sci') {
            const [mant, exp] = value.toExponential((spec.digits || 2) - 1).split('e');
            text = `${mant}\\times10^{${parseInt(exp, 10)}}`;
        }
        else if (spec.format === 'integer') text = String(value);
        else if (spec.format === 'pct') text = (value * 100).toPrecision(spec.digits || 3);
        else if (spec.format === 'sci') {
            const [mant, exp] = value.toExponential((spec.digits || 2) - 1).split('e');
            text = `${mant}&times;10<sup>${parseInt(exp, 10)}</sup>`;
        }
        else text = value.toPrecision(spec.digits || 3);
        return inMath ? text : `<span data-claim="${id}">${text}</span>`;
    });
}
module.exports = { renderClaims };
