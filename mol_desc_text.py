from openclatura import describe_human
import re
from collections import defaultdict

def describe_molecule(smiles):
    d = describe_human(smiles)
    text = d.text

    # Normalize whitespace
    text = re.sub(r"\s+", " ", text).strip()

    # Split into sentences
    sentences = re.split(r'(?<=[.])\s+', text)

    data = defaultdict(list)

    for s in sentences:
        if s.startswith("The molecule is named"):
            data["iupac_name"] = s.removeprefix("The molecule is named ").rstrip(".")

        elif s.startswith("The molecule is built around"):
            data["parent"] = s.removeprefix("The molecule is built around ").rstrip(".")

        elif s.startswith("Within that parent framework"):
            data["parent_sub_details"].append(
                s.removeprefix("Within that parent framework, ").rstrip(".")
            )

        elif s.startswith("The principal characteristic feature is"):
            data["principal_feature"] = (
                s.removeprefix("The principal characteristic feature is ")
                .rstrip(".")
            )

    try:
        parent = data['parent']
        principal_feature = data['principal_feature']
        return{            
            'parent': parent if parent else "No preferred Parent to show",
            'principal_feature': principal_feature if principal_feature else "No preferred Principal feature to show"
        }

        # print(data["parent"])
        # print(data["principal_feature"])
        # print(data["parent_sub_details"])
    
    except:
        return{            
            'parent': "No preferred Parent to show",
            'principal_feature': "No preferred Principal feature to show"
            }





# def describe_molecule(smiles):
#     d = describe_human(smiles)
#     text = d.text

#     data = {}

#     patterns = {
#         "iupac_name": r"The molecule is named (.+?)\.",
#         "parent": r"The molecule is built around (.+?)\.",
#         "principal_feature": r"The principal characteristic feature is (.+?)\.",
#     }

#     for key, pattern in patterns.items():
#         m = re.search(pattern, text)
#         if m:
#             data[key] = m.group(1)

#     try:
#         iupac_name = data['iupac_name'] 
#         parent = data['parent']
#         principal_feature = data['principal_feature']
#         return{
#             'iupac_name': iupac_name if iupac_name else "No preferred IUPAC name",
#             'parent': parent if parent else "No preferred Parent to show",
#             'principal_feature': principal_feature if principal_feature else "No preferred Principal feature to show"
#         }
    
#     except:
#         return{
#             'iupac_name': "No preferred IUPAC name",
#             'parent': "No preferred Parent to show",
#             'principal_feature': "No preferred Principal feature to show"
#             }
