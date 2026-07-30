# Release checklist for version 2.2.0

- [x] Preserve the verified author order, affiliations and ORCIDs.
- [x] Preserve MIT for code and CC BY 4.0 for original derived outputs,
  figures and documentation.
- [x] Exclude original GEO files, third-party supplements and author-provided
  raw objects.
- [x] Exclude the manuscript, cover letter and submission-only DOCX/PDF files.
- [x] Include all 23 frozen R scripts.
- [x] Include frozen figures, source data, supplementary results, audits and
  validation summaries.
- [x] Complete privacy, path, source-copy and syntax validation.
- [x] Generate and verify the public-release SHA-256 manifest.
- [x] Commit the reviewed public release candidate:
  `7e257ec4621e35d053dfed9fe1115c22ab2eb581`.
- [x] Push the release branch and merge PR #1 into `main`:
  `fda9b0a8d053399f08a09855f54e55c8a97005bd`.
- [x] Merge the release-state cleanup in PR #2:
  `d4c8af5a5f8d0c11cbefd062df57daa921aa3631`.
- [ ] Create GitHub tag and Release `v2.2.0`.
- [ ] Build the Zenodo archive from the exact tagged Git state.
- [x] Create a new version under Zenodo concept DOI
  `10.5281/zenodo.21677617` and reserve version DOI
  `10.5281/zenodo.21700808`.
- [ ] Verify the Zenodo file checksum, creators, licences, version and DOI.
- [ ] Publish the Zenodo version and test both DOI and GitHub URLs signed out.
- [ ] Only then update the manuscript and cover letter with the real version
  DOI and publication date.
- [ ] Rebuild and re-audit the final submission package.
