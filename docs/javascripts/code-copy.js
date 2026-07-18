document$.subscribe(() => {
  document.querySelectorAll("pre > code").forEach((code) => {
    const walker = document.createTreeWalker(code, NodeFilter.SHOW_TEXT);
    let lastText = null;
    while (walker.nextNode()) {
      lastText = walker.currentNode;
    }
    if (lastText?.data.endsWith("\n")) {
      lastText.data = lastText.data.slice(0, -1);
    }
  });
});
