# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Modified for genomic-variant-research from skills/gnomad_database/scripts/get_gene_constraint.py in
# https://github.com/google-deepmind/science-skills at commit 370f8f33e5bb27dffe8289ea9bdec436686613ce:
#   - renamed
#   - pinned polite-http==0.1.3
#   - --output may be a bare filename
#   - query also returns lof_z and mis_z

"""Fetches gene constraint metrics from gnomAD.

This script retrieves constraint metrics (pLI, LOEUF, etc.) for a given gene
from the gnomAD database using its GraphQL API. It enforces a rate limit
of 10 queries per minute to respect the API's usage policy.
"""

# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "polite-http==0.1.3",
# ]
# ///

import argparse
import json
import os

from polite_http import http_client

# Respect 10 queries per minute requirement.
CLIENT = http_client.HttpClient("https://gnomad.broadinstitute.org", qps=0.1666)


def get_gene_constraint(
    gene_symbol: str,
    output_path: str,
):
  """Fetch gene constraint from gnomAD."""
  url = "https://gnomad.broadinstitute.org/api"
  query = """
  query($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
    gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
      gene_id
      symbol
      gnomad_constraint {
        pli
        oe_lof
        oe_lof_lower
        oe_lof_upper
        oe_mis
        oe_mis_lower
        oe_mis_upper
        lof_z
        mis_z
      }
    }
  }
  """

  variables = {"geneSymbol": gene_symbol, "referenceGenome": "GRCh38"}

  response_data = CLIENT.fetch_json(
      url,
      method="POST",
      json_body={"query": query, "variables": variables},
  )
  result = json.dumps(response_data, indent=2)
  os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
  with open(output_path, "w") as f:
    f.write(result)
    f.write("\n")


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
      description="Fetch gene constraint from gnomAD"
  )
  parser.add_argument("--gene", required=True, help="Gene symbol (e.g. PCSK9)")
  parser.add_argument(
      "--output",
      "-o",
      required=True,
      help="Output file path. Prints to stdout if not specified.",
  )
  args = parser.parse_args()

  get_gene_constraint(args.gene, args.output)
