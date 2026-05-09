#!/usr/bin/env node

const fs = require('fs');
const path = require('path');

function readJsonIfExists(filePath) {
    if (!fs.existsSync(filePath)) {
        return null;
    }
    try {
        const raw = fs.readFileSync(filePath, 'utf8');
        try {
            return JSON.parse(raw);
        } catch (error) {
            const sanitized = raw
                .replace(/:\s*-Infinity\b/g, ': null')
                .replace(/:\s*Infinity\b/g, ': null')
                .replace(/:\s*NaN\b/g, ': null');
            return JSON.parse(sanitized);
        }
    } catch (error) {
        console.warn(`⚠️  Failed to parse JSON: ${filePath}`);
        return null;
    }
}

function formatFixedNumber(value, decimals) {
    if (value === null || value === undefined) {
        return '';
    }
    const num = Number(value);
    if (!Number.isFinite(num)) {
        return '';
    }
    return num.toFixed(decimals);
}

function formatSignedNumber(value, decimals) {
    const fixed = formatFixedNumber(value, decimals);
    if (!fixed) {
        return '';
    }
    const num = Number(value);
    return `${num >= 0 ? '+' : ''}${fixed}`;
}

function formatCI(lower, upper, decimals) {
    const lo = formatSignedNumber(lower, decimals);
    const hi = formatSignedNumber(upper, decimals);
    if (!lo || !hi) {
        return '';
    }
    return `[${lo}, ${hi}]`;
}

function formatIntegerWithCommas(value) {
    if (value === null || value === undefined) {
        return '';
    }
    const num = Number(value);
    if (!Number.isFinite(num)) {
        return '';
    }
    return Math.round(num).toLocaleString('en-US');
}

function formatIntegerLatex(value) {
    const s = formatIntegerWithCommas(value);
    if (!s) {
        return '';
    }
    return s.replace(/,/g, '{,}');
}

function formatPValueMantissaExp(pValue, mantissaDecimals = 1) {
    if (pValue === null || pValue === undefined) {
        return null;
    }
    const p = Number(pValue);
    if (!Number.isFinite(p) || p <= 0) {
        return null;
    }
    const exp = Math.floor(Math.log10(p));
    let mantissa = p / Math.pow(10, exp);
    let mantissaStr = mantissa.toFixed(mantissaDecimals);
    if (Number(mantissaStr) >= 10) {
        mantissa = mantissa / 10;
        mantissaStr = mantissa.toFixed(mantissaDecimals);
        return { mantissa: mantissaStr, exp: exp + 1 };
    }
    return { mantissa: mantissaStr, exp };
}

function formatPValueLatex(pValue, mantissaDecimals = 1) {
    const parts = formatPValueMantissaExp(pValue, mantissaDecimals);
    if (!parts) {
        return null;
    }
    return `${parts.mantissa} \\times 10^{${parts.exp}}`;
}

function formatPValueCell(pValue, options = {}) {
    const { includePrefix = false } = options;
    if (pValue === null || pValue === undefined) {
        return '';
    }
    const p = Number(pValue);
    if (!Number.isFinite(p)) {
        return '';
    }
    if (p >= 0.001) {
        const s = formatFixedNumber(p, 3);
        return includePrefix ? `p = ${s}` : s;
    }
    const latex = formatPValueLatex(p);
    if (!latex) {
        return '';
    }
    return includePrefix ? `p = ${latex}` : latex;
}

function formatPValueHtmlCell(pValue) {
    if (pValue === null || pValue === undefined) {
        return '';
    }
    const p = Number(pValue);
    if (!Number.isFinite(p) || p <= 0) {
        return '';
    }
    if (p >= 0.001) {
        return formatFixedNumber(p, 3);
    }
    const latex = formatPValueLatex(p);
    if (!latex) {
        return '';
    }
    return `$${latex}$`;
}

function safeGet(obj, pathParts) {
    let cur = obj;
    for (const part of pathParts) {
        if (cur === null || cur === undefined) {
            return undefined;
        }
        cur = cur[part];
    }
    return cur;
}

function createInjectionContext() {
    const resultsDir = path.join(__dirname, '..', 'results');

    const step005 = readJsonIfExists(path.join(resultsDir, 'step005_fitting_results.json'));

    if (!step005) {
        console.warn('⚠️  Missing step005_fitting_results.json; placeholders may remain unresolved.');
    }

    const ctx = {
        flyby: {
            beta: {},
            statistics: {},
            heterogeneity: {},
            bootstrap: {},
            loo: {}
        }
    };

    if (step005) {
        // Individual fits
        const fits = step005.individual_fits || {};
        for (const [key, fit] of Object.entries(fits)) {
            if (fit.fit && fit.fit.beta_fitted) {
                ctx.flyby.beta[key] = {
                    fitted: formatFixedNumber(fit.fit.beta_fitted, 10),
                    uncertainty: formatFixedNumber(fit.fit.beta_uncertainty, 10),
                    ppn_gamma: formatFixedNumber(fit.fit.ppn_gamma_deviation, 15),
                    ppn_valid: fit.fit.ppn_valid ? '✓' : '✗'
                };
            }
        }

        // Summary statistics
        const stats = step005.summary_statistics || {};
        const betaStats = stats.beta_statistics || {};
        ctx.flyby.statistics = {
            weighted_mean: formatFixedNumber(betaStats.weighted_mean, 10),
            weighted_uncertainty: formatFixedNumber(betaStats.weighted_uncertainty, 10),
            inflated_uncertainty: formatFixedNumber(betaStats.inflated_uncertainty, 10),
            mean: formatFixedNumber(betaStats.mean, 10),
            std: formatFixedNumber(betaStats.std, 10),
            recommended_beta: formatFixedNumber(stats.recommended_beta, 10),
            recommended_uncertainty: formatFixedNumber(stats.recommended_uncertainty, 10)
        };

        // Heterogeneity tests
        const hetero = stats.heterogeneity_tests || {};
        ctx.flyby.heterogeneity = {
            chi_squared: formatFixedNumber(hetero.chi_squared, 2),
            reduced_chi_squared: formatFixedNumber(hetero.reduced_chi_squared, 2),
            I_squared: formatFixedNumber(hetero.I_squared_percent, 3),
            cochran_Q: formatFixedNumber(hetero.cochran_Q, 2),
            interpretation: hetero.heterogeneity_interpretation || ''
        };

        // Bootstrap results
        const bootstrap = step005.robustness_analysis?.bootstrap || {};
        ctx.flyby.bootstrap = {
            n_bootstrap: formatIntegerWithCommas(bootstrap.n_bootstrap),
            mean: formatFixedNumber(bootstrap.bootstrap_mean, 10),
            median: formatFixedNumber(bootstrap.bootstrap_median, 10),
            std: formatFixedNumber(bootstrap.bootstrap_std, 10),
            ci_95_lower: formatFixedNumber(bootstrap.ci_95_lower, 10),
            ci_95_upper: formatFixedNumber(bootstrap.ci_95_upper, 10),
            ci_68_lower: formatFixedNumber(bootstrap.ci_68_lower, 10),
            ci_68_upper: formatFixedNumber(bootstrap.ci_68_upper, 10)
        };

        // Leave-one-out results
        const loo = step005.robustness_analysis?.leave_one_out || {};
        const looResults = loo.leave_one_out_results || {};
        for (const [key, result] of Object.entries(looResults)) {
            ctx.flyby.loo[key] = {
                beta_without: formatFixedNumber(result.beta_without_this, 10),
                remaining_n: formatIntegerWithCommas(result.remaining_n)
            };
        }
        ctx.flyby.loo.stability_coefficient = formatFixedNumber(loo.stability_coefficient, 2);
        ctx.flyby.loo.robust = loo.conclusion_robust ? 'Yes' : 'No';
        ctx.flyby.loo.interpretation = loo.interpretation || '';
    }

    return ctx;
}

function injectPlaceholders(template, context) {
    const unresolved = new Set();
    const replaced = template.replace(/\{\{\s*([^}]+?)\s*\}\}/g, (match, expr) => {
        const pathParts = expr.split('.').map(s => s.trim()).filter(Boolean);
        const value = safeGet(context, pathParts);
        if (value === undefined || value === null) {
            unresolved.add(expr);
            return match;
        }
        return String(value);
    });
    if (unresolved.size > 0) {
        console.warn(`⚠️  Unresolved placeholders (${unresolved.size}): ${Array.from(unresolved).slice(0, 10).join(', ')}${unresolved.size > 10 ? ', ...' : ''}`);
    }
    return replaced;
}

async function buildStaticSite() {
    console.log('🔨 Building static site...');
    
    try {
        // Clean dist directory
        const distDir = path.join(__dirname, 'dist');
        if (fs.existsSync(distDir)) {
            console.log('🧹 Cleaning dist directory...');
            fs.rmSync(distDir, { recursive: true, force: true });
        }
        
        // Read the manifest
        const manifestPath = path.join(__dirname, 'manifest.json');
        const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
        
        // Read the base index.html
        const indexPath = path.join(__dirname, 'index.html');
        let indexContent = fs.readFileSync(indexPath, 'utf8');
        
        const injectionContext = createInjectionContext();

        // Build the component content
        let componentsHtml = '';
        
        // Sort sections by order and load each component
        const sortedSections = manifest.sections.sort((a, b) => a.order - b.order);
        
        for (const section of sortedSections) {
            console.log(`📄 Loading section: ${section.title}`);
            
            const componentPath = path.join(__dirname, 'components', section.file);
            
            if (fs.existsSync(componentPath)) {
                const componentHtmlRaw = fs.readFileSync(componentPath, 'utf8');
                const componentHtml = injectPlaceholders(componentHtmlRaw, injectionContext);
                
                // Wrap component in section container (matching the dynamic loader)
                componentsHtml += `
                <section id="${section.id}" class="manuscript-section" data-section="${section.title}">
                    ${componentHtml}
                </section>`;
            } else {
                console.warn(`⚠️  Component not found: ${section.file}`);
                componentsHtml += `
                <section id="${section.id}" class="manuscript-section" data-section="${section.title}">
                    <div style="background-color: #ffe6e6; border: 1px solid #ff9999; padding: 15px; margin: 20px 0; border-radius: 5px;">
                        <h3 style="color: #cc0000; margin-top: 0;">Missing Section: ${section.title}</h3>
                        <p>Component file <code>components/${section.file}</code> not found</p>
                    </div>
                </section>`;
            }
        }
        
        // Replace the dynamic content with static content
        const staticContent = indexContent
            .replace(
                /<div id="loading".*?<\/div>\s*<div id="manuscript-content".*?<\/div>/s,
                () => `<div id="manuscript-content">${componentsHtml}</div>`
            )
            .replace(/<div id="loading"[^>]*>[\s\S]*?<\/div>/g, '')
            .replace(/<div id="manuscript-content"[^>]*style="display:\s*none;"[^>]*>[\s\S]*?<\/div>/g, '')
            .replace(
                /<!-- Component Loading Script -->[\s\S]*?<\/script>/g,
                `<!-- Static build - components pre-loaded -->\n    <script>\n        document.addEventListener("DOMContentLoaded", function() {\n            if (window.MathJax && window.MathJax.typesetPromise) {\n                window.MathJax.typesetPromise().catch(function (err) {\n                    console.error("MathJax error:", err.message);\n                });\n            }\n        });\n    </script>`
            )
            .replace(
                '<main id="main-content" role="main">',
                '<!-- This is a statically built version for SEO/deployment -->\n    <main id="main-content" role="main">'
            );
        
        // Create dist directory if it doesn't exist
        if (!fs.existsSync(distDir)) {
            fs.mkdirSync(distDir, { recursive: true });
        }
        
        // Write the built file
        const outputPath = path.join(distDir, 'index.html');
        fs.writeFileSync(outputPath, staticContent, 'utf8');
        
        // Copy necessary static assets to dist
        // EXPLICITLY COPY STYLES DIRECTORY
        const assetDirs = ['public', 'figures', 'data', 'styles'];
        for (const assetDir of assetDirs) {
            const srcPath = path.join(__dirname, assetDir);
            const destPath = path.join(distDir, assetDir);
            
            if (fs.existsSync(srcPath)) {
                console.log(`📁 Copying ${assetDir}/`);
                copyRecursiveSync(srcPath, destPath);
            }
        }

        // Copy simple-styles.css if it exists
        const simpleStylesSrc = path.join(__dirname, 'simple-styles.css');
        const simpleStylesDest = path.join(distDir, 'simple-styles.css');
        if (fs.existsSync(simpleStylesSrc)) {
            fs.copyFileSync(simpleStylesSrc, simpleStylesDest);
            console.log('📁 Copied simple-styles.css');
        }
        
        // Copy figures from results/figures/ and root-level results PNGs.
        const resultsFiguresPath = path.join(__dirname, '..', 'results', 'figures');
        const distFiguresPath = path.join(distDir, 'figures');
        const distPublicFiguresPath = path.join(distDir, 'public', 'figures');
        if (fs.existsSync(resultsFiguresPath)) {
            console.log('📁 Copying results/figures/ → dist/figures/');
            copyRecursiveSync(resultsFiguresPath, distFiguresPath);
            console.log('📁 Copying results/figures/ → dist/public/figures/');
            copyRecursiveSync(resultsFiguresPath, distPublicFiguresPath);
        }

        const resultsPath = path.join(__dirname, '..', 'results');
        if (fs.existsSync(resultsPath)) {
            const rootFigureFiles = fs.readdirSync(resultsPath)
                .filter((file) => /^step\d+_figure.*\.png$/i.test(file));

            for (const file of rootFigureFiles) {
                fs.mkdirSync(distFiguresPath, { recursive: true });
                fs.mkdirSync(distPublicFiguresPath, { recursive: true });
                fs.copyFileSync(
                    path.join(resultsPath, file),
                    path.join(distFiguresPath, file)
                );
                fs.copyFileSync(
                    path.join(resultsPath, file),
                    path.join(distPublicFiguresPath, file)
                );
            }

            if (rootFigureFiles.length > 0) {
                console.log(`📁 Copied ${rootFigureFiles.length} root results figure(s)`);
            }
        }
        
        // Copy manifest.json for reference
        fs.copyFileSync(manifestPath, path.join(distDir, 'manifest.json'));
        
        // Copy .nojekyll to dist root for GitHub Pages
        const nojekyllSrc = path.join(__dirname, 'public', '.nojekyll');
        const nojekyllDest = path.join(distDir, '.nojekyll');
        if (fs.existsSync(nojekyllSrc)) {
            fs.copyFileSync(nojekyllSrc, nojekyllDest);
            console.log('📁 Copied .nojekyll to dist root');
        } else {
            // Create it if it doesn't exist
            fs.writeFileSync(nojekyllDest, '');
            console.log('📁 Created .nojekyll in dist root');
        }

        // Copy robots.txt and sitemap.xml to dist root
        const rootFiles = ['404.html', 'robots.txt', 'sitemap.xml', 'CNAME', '29c6507763d2303d801cc8ed89d39f88.txt', 'favicon.ico'];
        for (const file of rootFiles) {
            const src = path.join(__dirname, 'public', file);
            const dest = path.join(distDir, file);
            if (fs.existsSync(src)) {
                fs.copyFileSync(src, dest);
                console.log(`📁 Copied ${file} to dist root`);
            }
        }
        
        // Copy citation files to dist root
        const citationFiles = ['CITATION.cff', 'CITATION.bib', 'citation.json', 'codemeta.json', 'README.md'];
        const optionalCitationFiles = new Set(['LICENSE']);
        const projectRoot = path.join(__dirname, '..');
        
        for (const file of [...citationFiles, ...optionalCitationFiles]) {
            const candidates = [
                path.join(projectRoot, file),
                path.join(__dirname, file),
            ];
            const src = candidates.find(candidate => fs.existsSync(candidate));
            const dest = path.join(distDir, file);
            if (src) {
                fs.copyFileSync(src, dest);
                console.log(`📁 Copied ${file} to dist root`);
            } else if (!optionalCitationFiles.has(file)) {
                console.warn(`⚠️  Missing citation file: ${file}`);
            }
        }
        
        // Copy .well-known directory
        const wellKnownSrc = path.join(__dirname, 'public', '.well-known');
        const wellKnownDest = path.join(distDir, '.well-known');
        if (fs.existsSync(wellKnownSrc)) {
            if (!fs.existsSync(wellKnownDest)) {
                fs.mkdirSync(wellKnownDest, { recursive: true });
            }
            const files = fs.readdirSync(wellKnownSrc);
            for (const file of files) {
                fs.copyFileSync(path.join(wellKnownSrc, file), path.join(wellKnownDest, file));
            }
            console.log(`📁 Copied .well-known/ to dist`);
        }
        
        // Markdown generation is now handled by dev-server.js to avoid duplication
        
        console.log('✅ Static site built successfully!');
        console.log(`📁 Output: ${outputPath}`);
        console.log(`📊 Generated ${manifest.sections.length} sections`);
        console.log('🚀 Ready for deployment');
        
    } catch (error) {
        console.error('❌ Build failed:', error);
        process.exit(1);
    }
}

// Helper function to copy directories recursively
function copyRecursiveSync(src, dest) {
    const exists = fs.existsSync(src);
    const stats = exists && fs.statSync(src);
    const isDirectory = exists && stats.isDirectory();
    
    if (isDirectory) {
        if (!fs.existsSync(dest)) {
            fs.mkdirSync(dest, { recursive: true });
        }
        fs.readdirSync(src).forEach(childItemName => {
            copyRecursiveSync(
                path.join(src, childItemName),
                path.join(dest, childItemName)
            );
        });
    } else {
        // Filter out raw data files
        if (src.endsWith('.csv') || src.endsWith('.dat') || src.endsWith('.nc')) {
            return;
        }
        fs.copyFileSync(src, dest);
    }
}

// Run if called directly
if (require.main === module) {
    buildStaticSite();
}

module.exports = { buildStaticSite };
