<h1 align="center">MetalFlare</h1>

<h4 align="center">Investigating metal-sensing green fluorescent protein</h4>

<p align="center">
    <a href="https://gitlab.com/oasci/metalflare/" target="_blank">
        <img alt="GitHub repo size" src="https://img.shields.io/github/repo-size/oasci/metalflare">
    </a>
    <a href="https://github.com/psf/black" target="_blank">
        <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Black style">
    </a>
</p>

TODO:

## Communication

> No research should be done alone.

We use [this repository's issues](https://gitlab.com/oasci/metalflare/issues) as our todo list.
Asynchronous conversations about any of the tasks should be included as issue comments.
Synchronous meetings happening in-person or virtually should have meeting minutes stored in the [appropriate directory](docs/01-management/03-meetings).

## Deploying

A note to maintainers.

We use [bump-my-version](https://github.com/callowayproject/bump-my-version) to release a new version.
This will create a git tag that is used by [poetry-dynamic-version](https://github.com/mtkennerly/poetry-dynamic-versioning) to generate version strings.

However, we are using [Calendar Versioning](https://calver.org/) which means we need to manually specify new versions.
For example, to bump the version to November 8, 2024, you would run the following command after activating the relevant conda environment.

```bash
bump-my-version bump --new-version 2024.11.8
```

After releasing a new version, you need to push and include all tags.

```bash
git push --follow-tags
```

## License

Code contained in this project is released under the [MIT License](https://spdx.org/licenses/MIT.html) as specified in [`LICENSE_CODE`](https://gitlab.com/oasci/metalflare/blob/main/LICENSE_CODE.md), which grants you the freedom to use, modify, and distribute it as long as you include the original copyright notice and disclaimer.

> Portions of this code were incorporated and adapted with permission from [metalflare](https://gitlab.com/oasci/metalflare) by OASCI under the [MIT License](https://gitlab.com/oasci/metalflare/blob/main/LICENSE_CODE.md).

All other data, information, documentation, and associated content provided within this project are released under the [Creative Commons Attribution-NoDerivatives 4.0 International License (CC BY-ND 4.0)](https://creativecommons.org/licenses/by-nd/4.0/) as specified in [`LICENSE_INFO`](https://gitlab.com/oasci/metalflare/blob/main/LICENSE_INFO.md).
You can freely share and use the material for any purpose, including commercial use, as long as you follow the license terms, which require attribution and prohibit distribution of modified material or imposing additional restrictions.

### Open-source release

On **January 1, 2033**, the above CC BY-ND 4.0 license is voided, removed, and superseded by [`LICENSE_INFO_OPEN`](https://gitlab.com/oasci/metalflare/blob/main/LICENSE_INFO_OPEN.md), and the following content.
This aforementioned date is subject to change to a time earlier than **January 1, 2033**, but never later.

> All other data, information, documentation, and associated content provided within this project are released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) as specified in [`LICENSE_INFO_OPEN`](https://gitlab.com/oasci/metalflare/blob/main/LICENSE_INFO_OPEN.md).
> This means you are free to share and adapt the non-code elements, but you must give appropriate credit to the original source and indicate if changes were made.
