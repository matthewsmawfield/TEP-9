"""Map pipeline flyby keys to IRI trajectory profile names."""

IRI_MISSION_MAP = {
    "NEAR": "NEAR_1998",
    "Galileo_1990": "Galileo_1990",
    "Galileo_1992": "Galileo_1992",
    "Cassini": "Cassini_1999",
    "Rosetta_2005": "Rosetta_2005",
    "Rosetta_2007": "Rosetta_2007",
    "Rosetta_2009": "Rosetta_2009",
    "MESSENGER_2005": "MESSENGER_2005",
    "Juno": "Juno_2013",
    "Stardust": "Stardust_2001",
}


def resolve_iri_mission(mission_name: str) -> str:
    return IRI_MISSION_MAP.get(mission_name, mission_name)
