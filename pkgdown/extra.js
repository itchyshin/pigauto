(function () {
  function enhancePage() {
    // Give downlit's one-character base-R c() links a useful accessible name
    // while keeping the visible signature and destination unchanged.
    document.querySelectorAll('main a[href="https://rdrr.io/r/base/c.html"]').forEach(function (link) {
      if (link.textContent.replace(/\s+/g, " ").trim() !== "c") return;
      if (!link.hasAttribute("aria-label")) {
        link.setAttribute("aria-label", "c() documentation (base R)");
      }
    });

    // The pkgdown skip link points to #main. Give that landmark a programmatic
    // tab stop so fragment navigation moves keyboard focus as well as scroll.
    var main = document.querySelector("main#main");
    if (main && !main.hasAttribute("tabindex")) main.setAttribute("tabindex", "-1");

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

    document.querySelectorAll("a.anchor[href^='#']").forEach(function (anchor) {
      var target = anchor.parentElement;
      if (!target || !target.matches("h1, h2, h3, h4, h5, h6, dt")) return;

      var labelSource = target.cloneNode(true);
      labelSource.querySelectorAll("a.anchor").forEach(function (link) {
        link.remove();
      });
      var label = labelSource.textContent.replace(/\s+/g, " ").trim();
      if (label) anchor.setAttribute("aria-label", "Link to " + label);
    });

    // Keep the type comparison tables in predict.pigauto_fit readable on
    // narrow screens. The roxygen HTML tables are intentionally semantic
    // tables; add only the responsive container used by the shared cue code.
    ["pooling-decoded-draws", "what-se-means-per-trait-type"].forEach(function (id) {
      var heading = document.getElementById(id);
      var section = heading && heading.closest(".section");
      var table = section && section.querySelector("table");
      if (!table || table.parentElement.classList.contains("table-responsive")) return;
      var wrapper = document.createElement("div");
      wrapper.className = "table-responsive";
      table.parentNode.insertBefore(wrapper, table);
      wrapper.appendChild(table);
    });

    var cueNumber = 0;
    function updateScrollCues() {
      document.querySelectorAll("main pre, main .table-responsive, main .math-overflow-region").forEach(function (target) {
        var isTable = target.classList.contains("table-responsive");
        var isEquation = target.classList.contains("math-overflow-region");
        var overflows = target.scrollWidth > target.clientWidth + 2;
        var cue = target.previousElementSibling;
        if (overflows && !(cue && cue.classList.contains("pa-scroll-cue"))) {
          cueNumber += 1;
          cue = document.createElement("p");
          cue.className = "pa-scroll-cue";
          cue.id = "pa-scroll-cue-" + cueNumber;
          cue.textContent = isTable
            ? "This table scrolls horizontally to show all columns."
            : (isEquation
              ? "This equation scrolls horizontally to show all terms."
              : "This code block scrolls horizontally to show the full line.");
          target.parentNode.insertBefore(cue, target);
          target.classList.add("pa-has-overflow");
          if (!target.hasAttribute("tabindex")) {
            target.setAttribute("tabindex", "0");
            target.setAttribute("data-pa-added-tabindex", "true");
          }
          var describedBy = (target.getAttribute("aria-describedby") || "").split(/\s+/).filter(Boolean);
          if (describedBy.indexOf(cue.id) === -1) describedBy.push(cue.id);
          target.setAttribute("aria-describedby", describedBy.join(" "));
        } else if (!overflows && cue && cue.classList.contains("pa-scroll-cue")) {
          var cueId = cue.id;
          cue.remove();
          target.classList.remove("pa-has-overflow");
          if (target.getAttribute("data-pa-added-tabindex") === "true") {
            target.removeAttribute("tabindex");
            target.removeAttribute("data-pa-added-tabindex");
          }
          var remaining = (target.getAttribute("aria-describedby") || "").split(/\s+/).filter(function (id) {
            return id && id !== cueId;
          });
          if (remaining.length) target.setAttribute("aria-describedby", remaining.join(" "));
          else target.removeAttribute("aria-describedby");
        }
      });
    }
    updateScrollCues();
    window.addEventListener("load", updateScrollCues, { once: true });
    if (document.fonts && document.fonts.ready) {
      document.fonts.ready.then(updateScrollCues);
    }
    if (window.ResizeObserver) {
      var resizeObserver = new ResizeObserver(updateScrollCues);
      resizeObserver.observe(document.querySelector("main") || document.body);
      document.querySelectorAll("main .math-overflow-region").forEach(function (target) {
        resizeObserver.observe(target);
      });
    } else {
      window.addEventListener("resize", updateScrollCues);
    }
    if (window.MutationObserver) {
      document.querySelectorAll("main .math-overflow-region").forEach(function (target) {
        new MutationObserver(updateScrollCues).observe(target, { childList: true, subtree: true, attributes: true });
      });
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", enhancePage, { once: true });
  } else {
    enhancePage();
  }
})();
