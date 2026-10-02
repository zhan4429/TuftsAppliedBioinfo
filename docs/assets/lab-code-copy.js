(function () {
  const copyIcon = `
    <svg viewBox="0 0 24 24" aria-hidden="true" focusable="false">
      <rect x="8" y="8" width="11" height="11" rx="2"></rect>
      <path d="M5 15H4a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h8a2 2 0 0 1 2 2v1"></path>
    </svg>
  `;
  const checkIcon = `
    <svg viewBox="0 0 24 24" aria-hidden="true" focusable="false">
      <path d="M20 6 9 17l-5-5"></path>
    </svg>
  `;

  function getCodeText(pre) {
    const code = pre.querySelector("code");
    return code ? code.textContent || "" : pre.textContent || "";
  }

  async function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      await navigator.clipboard.writeText(text);
      return;
    }

    const textarea = document.createElement("textarea");
    textarea.value = text;
    textarea.setAttribute("readonly", "");
    textarea.style.position = "fixed";
    textarea.style.inset = "0 auto auto 0";
    textarea.style.opacity = "0";
    document.body.append(textarea);
    textarea.select();
    document.execCommand("copy");
    textarea.remove();
  }

  function setButtonState(button, state) {
    if (state === "copied") {
      button.innerHTML = checkIcon;
      button.classList.add("is-copied");
      button.setAttribute("aria-label", "Copied");
      button.title = "Copied";
      return;
    }

    button.innerHTML = copyIcon;
    button.classList.remove("is-copied");
    button.setAttribute("aria-label", "Copy code");
    button.title = "Copy code";
  }

  function enhanceCodeBlock(pre) {
    if (pre.closest(".lab-code-block")) {
      return;
    }

    const wrapper = document.createElement("div");
    wrapper.className = "lab-code-block";
    pre.before(wrapper);
    wrapper.append(pre);

    const button = document.createElement("button");
    button.className = "lab-code-copy";
    button.type = "button";
    setButtonState(button);
    wrapper.append(button);

    button.addEventListener("click", async () => {
      try {
        await copyText(getCodeText(pre));
        setButtonState(button, "copied");
        window.setTimeout(() => setButtonState(button), 1400);
      } catch (_error) {
        button.setAttribute("aria-label", "Copy failed");
        button.title = "Copy failed";
      }
    });
  }

  document.addEventListener("DOMContentLoaded", () => {
    document.querySelectorAll(".lab-article pre.lab-code").forEach(enhanceCodeBlock);
  });
})();
