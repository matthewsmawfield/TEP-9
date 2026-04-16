const fs = require('fs');
const path = require('path');

const COMPONENTS_DIR = path.join(__dirname, '../site/components');
const OUTPUT_FILE = path.join(__dirname, '../15manuscript-tep-wb.md');

// Define the order of components
const components = [
    '1_introduction.html',
    '2_theory.html',
    '3_data.html',
    '4_results.html',
    '5_environment.html',
    '6_discussion.html',
    '7_conclusion.html'
];

let markdown = `# The Wide Binary Anomaly: Resolving the Gaia DR3 Controversy via Temporal Screening\n\n`;

components.forEach(file => {
    const filePath = path.join(COMPONENTS_DIR, file);
    if (fs.existsSync(filePath)) {
        let content = fs.readFileSync(filePath, 'utf8');
        
        // Simple HTML to Markdown conversion for headers and paragraphs
        content = content.replace(/<h2>(.*?)<\/h2>/g, '## $1');
        content = content.replace(/<h3>(.*?)<\/h3>/g, '### $1');
        content = content.replace(/<p>(.*?)<\/p>/g, '$1\n');
        content = content.replace(/<strong>(.*?)<\/strong>/g, '**$1**');
        content = content.replace(/<em>(.*?)<\/em>/g, '*$1*');
        content = content.replace(/<ol>/g, '');
        content = content.replace(/<\/ol>/g, '');
        content = content.replace(/<li>(.*?)<\/li>/g, '1. $1');
        
        markdown += content + '\n\n';
    } else {
        console.warn(`Warning: Component /Users/matthewsmawfield/www/TEP-JWST/zenodo.txt not found.`);
    }
});

fs.writeFileSync(OUTPUT_FILE, markdown);
console.log(`Generated ${OUTPUT_FILE}`);
