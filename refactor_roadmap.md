# pydna maintainability roadmap

Rule: legacy never breaks. Add a deprecation warning now, remove later.
Order: bloodless -> heaviest. One step = one PR.

```
  bloodless
     |
     v
  [PR1]  declare the public API with __all__      (no behavior change)
     |
     v
  [PR2]  make assembly2 the canonical Assembly     (old path still works, warns)
     |
     v
  [PR3]  one namespace to import the core types    (old imports still work)
     |
     v
  [PR4]  method registry: add a method = add a file
     |
     v
  [PR5]  split Dseq / Dseqrecord into mixins
     |
     v
  heaviest

  side track:  Makefile + CI share one command     (dev only, no runtime change)

  every PR:  tests pass, no public name removed without a warning
```

Now: PR1.
