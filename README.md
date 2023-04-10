## SHA-256-in-python

The SHA-256 coded in a straightforward way in python

#### Content
- The official [publication](http://csrc.nist.gov/publications/fips/fips180-4/fips-180-4.pdf) of SHA algorithm.
- A python numpy implementation of SHA-256 following the exact steps and notations as in the publication, in `sha256py`.

#### Tests
- `python3 sha256.py` will test the algorithm against hashlib 

#### Usage

```python
import sha256

message = list(b"some text")

digest=sha256.sha256(message)

hashedMessage=sha256.hexdigest(digest)
```
