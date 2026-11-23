# Step 7 - renumber references by order of first citation.
# Two references were appended at the tail, so the numbering has to be rebuilt
# from where each number first appears in the running text.

import io
import re

P = (r"F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/Manuscript_GutHeart_8genes"
     r"\Manuscript_GutHeart_8genes_v1.md")

t = io.open(P, encoding="utf-8").read()
body, sep, refs = t.partition("\n# References")
assert sep, "reference heading not found"

# ---- parse the existing entries ------------------------------------------
entries = {}
order_present = []
for ln in refs.split("\n"):
    m = re.match(r"^(\d+)\.\s+(.*\S)\s*$", ln)
    if m:
        n = int(m.group(1))
        entries[n] = m.group(2)
        order_present.append(n)

print("entries found: %d -> %s" % (len(entries), order_present))
assert order_present == sorted(order_present), "reference list is not ordered"

# ---- first-citation order in the body ------------------------------------
first = []
CITE = re.compile(r"\[(\d+(?:\s*,\s*\d+)*)\]")
for m in CITE.finditer(body):
    for k in m.group(1).split(","):
        n = int(k.strip())
        if n not in first:
            first.append(n)

print("first-citation order: %s" % first)
missing = [n for n in first if n not in entries]
extra = [n for n in entries if n not in first]
assert not missing, "cited but absent from the list: %s" % missing
assert not extra, "listed but never cited: %s" % extra

mapping = {old: i + 1 for i, old in enumerate(first)}
if all(k == v for k, v in mapping.items()):
    print("numbering already follows first citation; nothing to do")
else:
    print("remapping: %s" % mapping)

    def repl(m):
        nums = [int(x.strip()) for x in m.group(1).split(",")]
        return "[" + ",".join(str(mapping[n]) for n in nums) + "]"

    body = CITE.sub(repl, body)

    new_entries = [entries[old] for old in first]
    out = ["", "# References", ""]
    for i, txt in enumerate(new_entries, 1):
        out.append("%d. %s" % (i, txt))
    t = body + "\n" + "\n".join(out) + "\n"
    io.open(P, "w", encoding="utf-8", newline="\n").write(t)
    print("rewritten with %d references" % len(new_entries))

# ---- verification ---------------------------------------------------------
t = io.open(P, encoding="utf-8").read()
body, _, refs = t.partition("\n# References")
nums = [int(m.group(1)) for m in re.finditer(r"^(\d+)\.\s", refs, re.M)]
seen = []
for m in CITE.finditer(body):
    for k in m.group(1).split(","):
        n = int(k.strip())
        if n not in seen:
            seen.append(n)
print()
print("list numbers   :", nums)
print("first citations:", seen)
print("contiguous 1..n:", nums == list(range(1, len(nums) + 1)))
print("citations in ascending first-use order:", seen == sorted(seen))
print("every entry cited:", sorted(seen) == nums)
