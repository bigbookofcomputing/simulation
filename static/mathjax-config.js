window.MathJax = {
  tex: {
    inlineMath: [['$', '$'], ['\\(', '\\)']],
    displayMath: [['$$', '$$'], ['\\[', '\\]']],
    processEscapes: true,
    processEnvironments: true
  },

  options: {
    // Let MathJax process elements wrapped by arithmatex
    processHtmlClass: 'arithmatex',
    ignoreHtmlClass: 'tex2jax_ignore',

    // Don't skip DIVs (admonitions are DIVs)
    skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
  },

  // Fix: convert <span class="arithmatex"> back into inline math
  // so MathJax can render it (necessary for admonitions, tabs, details)
  renderActions: {
    find_arithmatex: [
      10,
      function () {
        document
          .querySelectorAll('span.arithmatex')
          .forEach(function (node) {
            const content = node.textContent;

            // Preserve original TeX if already in \( \) or \[ \]
            const math =
              content.startsWith('\\(') ||
              content.startsWith('\\[') ||
              content.startsWith('$$')
                ? content
                : '$' + content + '$';

            const text = document.createTextNode(math);
            node.replaceWith(text);
          });
      },
      ''
    ]
  }
};

// Re-run MathJax on each page load (required by mkdocs-material)
document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
