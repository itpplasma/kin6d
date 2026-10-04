# External literature corpus

The literature corpus is a separate reproducible research asset. Papers,
attachments, Zotero exports, DOI/download manifests, discovery queues and
metadata outputs live in Zotero and an external local workspace. kin6d stores
benchmark source identifiers and interpretation, not PDF papers or a corpus
dump. The [atlas](benchmarks.md) links formulas to the evidence they need.

Seed the external workspace with the supplied first-ring DOI list, a manual
edition/no-verified-DOI queue, and curated second-ring benchmark/code sources.
Resolve DOI metadata before interpreting descriptive titles as exact
bibliographic titles. Corrections to DOI/title pairs require recorded review;
never silently replace one work with another. Older books/translations can
have edition-specific identifiers even when first placed in the manual queue.

The idempotent workflow resolves Crossref/OpenAlex metadata, records legal OA
locations, downloads only public OA/repository candidates, validates the PDF,
and hashes accepted bytes. It uses no publisher login, cookies, authenticated
PDF scraping, or paywall bypass. Keep the downloaded URL, OA evidence,
retrieval time, SHA-256, metadata provenance, and review state. Failures and
unresolved identity or edition questions appear in `manual.tsv`.

Emit `downloaded.tsv`, `manual.tsv`, `metadata.jsonl`, and a Zotero-import
batch outside the solver repository. Check existing Zotero records before
creation; use the supported Web API for item/collection/attachment changes,
never mutate its SQLite database. Preserve user annotations and existing
collections/tags. Attachment promotion records its parent key and hash so
rerunning does not upload it again. Missing credentials are an explicit
blocked upload step, not a reason to fabricate a successful import.

References and forward citations feed a ranked review queue. Record the
root/relation, identifier, title, taxonomy match, metadata evidence and review
decision. Fetch bounded pages and record truncation/retry state. Discovery
does not admit or download papers automatically. Curate new roots before
running the resolver on them.

Collections cover foundations, exact benchmarks, collision operators,
relativistic/BBGKY, GC/DK/GK, neoclassical theory, numerical methods, structure,
certified numerics, codes/performance, real devices, and MHD/hybrid models.
Orthogonal tags describe equations, operators, geometries, species models,
proof methods, devices and codes. Collections organize reading; tags provide
cross-cutting retrieval. A paper's availability, peer review and relevance
do not prove a benchmark theorem; record its mathematical applicability in
the atlas before promoting the case.

Chris&AI
