import Util.user_input as UI
import pdb

ui = UI.client()

while True:
    key = ui.get_input()
    if key is not None:
        print key
