<h1 align="center">MetalFlare</h1>

<h4 align="center">Investigating metal-sensing green fluorescent protein</h4>

TODO:

## Communication

> No research should be done alone.

We use [this repository's issues](https://github.com/oasci/metalflare/issues) as our todo list.
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

By accessing or using the source code herein, you agree to the terms and conditions set forth in [`LICENSE_CODE`](https://github.com/oasci/metalflare/blob/main/LICENSE_CODE.md).
While you can use and distribute the provided source code, certain restrictions apply.
No warranties are provided, and neither party is liable for consequential damages.
Breach of terms may result in license termination, and both parties agree to indemnify each other against specific damages.
Please review [`LICENSE_CODE`](https://github.com/oasci/metalflare/blob/main/LICENSE_CODE.md) for comprehensive details.
This summary is provided for convenience, but the full license is the legally binding agreement that outlines the specifics of usage, restrictions, warranties, liabilities, and other important details.

All other data, information, documentation, and associated content provided within this project are released under the [Creative Commons Attribution-NoDerivatives 4.0 International License (CC BY-ND 4.0)](https://creativecommons.org/licenses/by-nd/4.0/) as specified in [`LICENSE_INFO`](https://github.com/oasci/metalflare/blob/main/LICENSE_INFO.md).
You can freely share and use the material for any purpose, including commercial use, as long as you follow the license terms, which require attribution and prohibit distribution of modified material or imposing additional restrictions.

These dual, semi-restrictive licenses ensure a balance between protecting our ongoing work and promoting open science while encouraging collaboration and proper attribution.

### Open-source release

On **January 1, 2033**, the above licenses and paragraphs are voided, removed, and superseded by [`LICENSE_CODE_OPEN`](https://github.com/oasci/metalflare/blob/main/LICENSE_CODE_OPEN.md), [`LICENSE_INFO_OPEN`](https://github.com/oasci/metalflare/blob/main/LICENSE_INFO_OPEN.md), and the following content.
This aforementioned date is subject to change to a time earlier than **January 1, 2033**, but never later.

> Code contained in this project is released under the [MIT License](https://spdx.org/licenses/MIT.html) as specified in [`LICENSE_CODE_OPEN`](https://github.com/oasci/metalflare/blob/main/LICENSE_CODE_OPEN.md), which grants you the freedom to use, modify, and distribute it as long as you include the original copyright notice and disclaimer.
>
> All other data, information, documentation, and associated content provided within this project are released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) as specified in [`LICENSE_INFO_OPEN`](https://github.com/oasci/metalflare/blob/main/LICENSE_INFO_OPEN.md).
> This means you are free to share and adapt the non-code elements, but you must give appropriate credit to the original source and indicate if changes were made.
>
> These dual licenses ensure a balance between open-source software and data accessibility while encouraging collaboration and proper attribution.
