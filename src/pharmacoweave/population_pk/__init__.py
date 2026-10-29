"""Population PK is not implemented.

Use an external NLME engine later, out of process. Importing this module does not run an analysis.
"""

from pharmacoweave.core.roadmap import status

STATUS = status("population_pk")
