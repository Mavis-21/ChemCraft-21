import streamlit as st
import sqlite3
import requests
import sys
import os

# --- NUCLEAR OPTION: DOWNLOAD LIBRARY DIRECTLY ---
# If streamlit-option-menu is missing, we download the raw code file 
# and load it manually. This bypasses pip/installation entirely.
try:
    from streamlit_option_menu import option_menu
except ImportError:
    print("⚠️ option_menu not found. downloading manually...")
    url = "https://raw.githubusercontent.com/victoryhb/streamlit-option-menu/master/streamlit_option_menu/__init__.py"
    r = requests.get(url)
    with open("streamlit_option_menu.py", "w") as f:
        f.write(r.text)
    # Now import from the file we just downloaded
    import streamlit_option_menu
    from streamlit_option_menu import option_menu

# ------------------------------- REST OF YOUR APP ---------------------------------

import py3Dmol
import time as t
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

# -------------------------------- FUNCTIONS ---------------------------------------

def initialisation():
    if "logged_in" not in st.session_state:
        st.session_state.logged_in = False
        st.session_state.guest = False
        st.session_state.user = None
    if "guest_history" not in st.session_state:
        st.session_state.guest_history = []
        st.session_state.count = 0

def update_page(x):
    if x == 0:
        st.session_state.page = "home"
    elif x == 1:
        st.session_state.page = "login"
    elif x == 2:
        st.session_state.page = "signup"

def update_accStatus(x):
    if x == 0:
        st.session_state.logged_in = True
    elif x == 1:
        st.session_state.guest = True

def update_mainpage(x):
    if x == 0:
        st.session_state.mainpage = "new"
    elif x == 1:
        st.session_state.mainpage = "aboutus"
    elif x == "sidebar":
        st.session_state.mainpage = "history"
    elif x == "sidebar_admin":
        st.session_state.mainpage = "Admin_tables"
    elif x == 2:
        st.session_state.mainpage = "Admin"

def toggle_fhistory(x):
    if x == 0:
        st.session_state.fhistory = False
    elif x == 1:
        st.session_state.fhistory = True

def toggle_afhistory(x):
    if x == 0:
        st.session_state.afhistory = False
    elif x == 1:
        st.session_state.afhistory = True

# --- SQLITE CONNECTION (ZERO CONFIG DATABASE) ---
def connection():
    con = sqlite3.connect('chemcraft.db', check_same_thread=False)
    cur = con.cursor()
    return con, cur

def users():
    con, cur = connection()
    try:
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='users'")
        if not cur.fetchone():
            return []
        cur.execute("select username from users")
        l = [i[0] for i in cur.fetchall()]
        return l
    except:
        return []
    finally:
        con.close()

def passwd_checker(u, p):
    con, cur = connection()
    try:
        cur.execute("select passwd from users where username=?", (u,))
        result = cur.fetchone()
        if result and result[0] == p:
            return True
        return False
    except:
        return False
    finally:
        con.close()

def user_table_exists():
    if st.session_state.user:
        con, cur = connection()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (st.session_state.user,))
        exists = cur.fetchone()
        con.close()
        return True if exists else False
    return False

def is_admin():
    if 'user' in st.session_state and st.session_state.user:
        con, cur = connection()
        try:
            cur.execute("SELECT typ from users where username = ?", (st.session_state.user,))
            r = cur.fetchone()
            if r and r[0] == "Admin":
                return True
        except:
            pass
        finally:
            con.close()
        return False
    return None

# -------------------------------------- SQL -------------------------------------------
def create_tables():
    con, cur = connection()
    q1 = """
            CREATE TABLE IF NOT EXISTS users(
                userid INTEGER PRIMARY KEY AUTOINCREMENT,
                username TEXT UNIQUE,
                passwd TEXT,
                email TEXT,
                typ TEXT
            )
            """
    cur.execute(q1)
    con.commit()
    con.close()

def create_usertable(user):
    con, cur = connection()
    # Using sanitized username as table name
    q = f"""create table if not exists "{user}" (
       user text,
       searched text,
       smiles text
       )""" 
    cur.execute(q)
    con.commit()
    con.close()

def get_history(username, x):
    if not username: return []
    con, cur = connection()
    try:
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (username,))
        if not cur.fetchone(): return []

        col = "searched" if x == "searched" else "smiles"
        cur.execute(f'SELECT {col} FROM "{username}"')
        result = cur.fetchall()
        return [i[0] for i in result] if result else []
    except:
        return []
    finally:
        con.close()

def get_tables():
    con, cur = connection()
    cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
    r = cur.fetchall()
    con.close()
    return [i[0] for i in r]

def update_history(iupac_input):
    smiles = iupac_to_smiles(iupac_input)
    if st.session_state.get("user") and smiles:
        current = get_history(st.session_state.user, "smiles")
        if smiles not in current:
            con, cur = connection()
            create_usertable(st.session_state.user)
            query = f'INSERT into "{st.session_state.user}" values(?, ?, ?)'
            cur.execute(query, (st.session_state.user, iupac_input, smiles))
            con.commit()
            con.close()

# _____________________________________ PAGES ________________________________________

def login():
    st.title("Welcome back to Chemcraft")
    st.header("Login")
    with st.form(key='login'):
        user = st.text_input("Username:", placeholder="Username")
        passwd = st.text_input("Password:", placeholder="Password", type='password')
        if st.form_submit_button("Login"):
            if not user or not passwd:
                st.warning("Please enter both username and password")
            elif user not in users():
                st.error("User does not exist")
            elif not passwd_checker(user, passwd):
                st.error("Password is incorrect")
            else:
                st.success("Logged in successfully")
                update_accStatus(0)
                st.session_state.user = user
                st.session_state.page = "dashboard"
                start() # Rerun
    st.button("Home", key='loginbutton', on_click=update_page, args=(0,))

def sign_up():
    st.title("Welcome to Chemcraft")
    st.header("Sign Up")
    with st.form(key='sign up'):
        user = st.text_input("Username:", placeholder="Username")
        passwd = st.text_input("Password:", placeholder="Password", type='password')
        email = st.text_input("Email:", placeholder="Email")
        typ = st.selectbox("Who are you?", ("High school Student", "College Student", "Professor", "Enthusiast", "Admin"))

        if st.form_submit_button("Sign Up"):
            if not user or not passwd:
                st.warning("Please fill in all fields")
            elif user in users():
                st.error("Username already exists")
            else:
                con, cur = connection()
                cur.execute("INSERT INTO users (username, passwd, email, typ) VALUES (?, ?, ?, ?)", (user, passwd, email, typ))
                con.commit()
                con.close()
                st.success("Account created")
                update_accStatus(0)
                st.session_state.user = user
                st.session_state.page = "dashboard"
                start()
    st.button("Home", key='signup', on_click=update_page, args=(0,))

def home():
    # STYLISH UI
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Orbitron:wght@800&display=swap');
    .stApp { background: radial-gradient(circle at 18% 23%, #121e30 0, #182033 60%, #141d28 100%); }
    .chem-header {
        font-family: 'Orbitron', sans-serif; font-size: 3.6rem; font-weight: 900; text-align: center;
        background: linear-gradient(110deg,#27e3fd 25%, #a855f7 60%, #ffb300 85%);
        -webkit-background-clip: text; -webkit-text-fill-color: transparent; margin: 36px 0 68px 0;
    }
    .entry-block {
        display: flex; gap: 0.9rem; margin: 43px 0; padding: 0 1rem;
        border-left: 6px solid #21e7abaf; border-radius: 1.3em 0 0 2.2em;
    }
    .entry-block-alt { border-left: 6px solid #4e9cffbf; margin-left: 85px; }
    .entry-bio { font-size: 1.18rem; color: #53fcf1; font-weight: bold; }
    .entry-intro { color: #f0f6ff; font-size: 1.03rem; }
    </style>
    """, unsafe_allow_html=True)

    st.markdown('<div class="chem-header">ChemCraft</div>', unsafe_allow_html=True)

    st.markdown('<div class="entry-block"><div><div class="entry-bio">🌟 New user?</div><div class="entry-intro">Create a free ChemCraft account.</div></div></div>', unsafe_allow_html=True)
    if st.button("Sign Up", key="signup_btn", type="primary"): update_page(2)

    st.markdown('<div class="entry-block entry-block-alt"><div><div class="entry-bio">🔐 Member?</div><div class="entry-intro">Log in to continue.</div></div></div>', unsafe_allow_html=True)
    if st.button("Log In", key="login_btn", type="secondary"): update_page(1)

    st.markdown('<div class="entry-block"><div><div class="entry-bio">👀 Guest?</div><div class="entry-intro">Try instantly.</div></div></div>', unsafe_allow_html=True)
    if st.button("Guest Mode", key="guest_btn"): update_accStatus(1)

# ___________________________________ RENDERING ___________________________________

def iupac_to_smiles(iupac_name):
    try:
        url = f"https://opsin.ch.cam.ac.uk/opsin/{iupac_name}.json"
        r = requests.get(url, timeout=5)
        if r.status_code == 200: return r.json().get("smiles", None)
    except:
        return None
    return None

def fetch_3d_structure(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToMolBlock(mol)
    except:
        return None

def rendering(iupac_input):
    if not iupac_input: return
    with st.spinner("Converting..."):
        smiles = iupac_to_smiles(iupac_input)
    
    if smiles:
        st.success(f"✅ SMILES: {smiles}")
        sdf_data = fetch_3d_structure(smiles)
        if sdf_data:
            # Simple 3D Viewer using Py3DMol
            view = py3Dmol.view(width=700, height=500)
            view.addModel(sdf_data, "mol")
            view.setStyle({'stick': {}})
            view.zoomTo()
            # Generate HTML for viewer
            html = view._make_html()
            st.components.v1.html(html, height=500)
        else:
            st.error("Could not generate structure.")
    else:
        st.error("Invalid IUPAC name.")

def fullrendering():
    st.title("🧪 3D Molecule Viewer")
    iupac_input = st.text_input("Enter IUPAC Name:")
    update_history(iupac_input)
    rendering(iupac_input)

# ___________________________________ SIDEBAR & MAIN _______________________________

def page_main():
    user = st.session_state.get("user")
    typ = "Admin" if is_admin() else "Civillian"
    
    with st.sidebar:
        st.button("New chat", on_click=update_mainpage, args=(0,))
        st.button("About Us", on_click=update_mainpage, args=(1,))
        
        # Try to use option_menu (which we force-downloaded if needed)
        try:
            # Simple History
            h = get_history(user, "searched") if st.session_state.logged_in else st.session_state.get("guest_history", [])
            if h:
                selected = option_menu("History", h[-5:], menu_icon="clock", default_index=0)
                if selected:
                    st.session_state.mainpage = "history"
                    st.session_state.sidebar = selected
        except Exception as e:
            st.error(f"Menu Error: {e}")

    if st.session_state.get("mainpage") == "new":
        fullrendering()
    elif st.session_state.get("mainpage") == "history":
        rendering(st.session_state.get("sidebar"))
    elif st.session_state.get("mainpage") == "aboutus":
        st.title("About ChemCraft")
        st.write("A chemistry visualization tool.")
    elif st.session_state.get("mainpage") == "Admin_tables":
        display_table(st.session_state.get("sidebar_admin"),0)

def start():
    if not(st.session_state.logged_in or st.session_state.guest):
        if "page" not in st.session_state: st.session_state.page = "home"
        
        if st.session_state.page == "home": home()
        elif st.session_state.page == "signup": sign_up()
        elif st.session_state.page == "login": login()
    else:
        st.rerun()

# ___________________________________ APP ENTRY ___________________________________

create_tables()
initialisation()

def Main():
    if st.session_state.logged_in or st.session_state.guest:
        page_main()
    else:
        start()

if __name__ == "__main__":
    Main()
