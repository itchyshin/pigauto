(function () {
  document.querySelectorAll("blockquote p").forEach(function (p) {
    var html = p.innerHTML;
    var m = html.match(/^\[!(WARNING|NOTE|IMPORTANT|TIP|CAUTION)\]\s*/i);
    if (!m) return;
    p.innerHTML = html.replace(/^\[![A-Z]+\]\s*/i, "");
    var bq = p.closest("blockquote");
    if (bq) {
      bq.classList.add("pa-callout", "pa-callout-" + m[1].toLowerCase());
    }
  });
})();
