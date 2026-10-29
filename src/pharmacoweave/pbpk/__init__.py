"""PBPK is not implemented.

Use an external engine adapter later, out of process. Importing this module does not run an analysis.
"""

from pharmacoweave.core.roadmap import status

STATUS = status("pbpk")
