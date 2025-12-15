import xml.etree.ElementTree as ET

# XMLファイル読み込み
tree = ET.parse("desc2025.xml")
root = tree.getroot()

# B-LymphocytesやHemic System / Blood CellsのDescriptorNameを検索
target_terms = ["B-Lymphocytes", "Blood Cells", "Hemic System"]

def get_mesh_ids(node, found_ids=set()):
    # DescriptorName
    name = node.findtext("DescriptorName/String")
    # DescriptorUI
    did = node.findtext("DescriptorUI")
    # TreeNumberListの取得
    tree_numbers = [tn.text for tn in node.findall("TreeNumberList/TreeNumber")]
    if any(tn.startswith("D12") for tn in tree_numbers):  # D12はHemic System系
        found_ids.add(did)
    # B-LymphocytesはD001711
    if name in target_terms:
        found_ids.add(did)
    return found_ids

mesh_ids = set()
for descriptor in root.findall("DescriptorRecord"):
    mesh_ids = get_mesh_ids(descriptor, mesh_ids)

print(mesh_ids)

